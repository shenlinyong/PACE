# PACE 模型与公式

[英文公式](FORMULA.md) · [数据准备](PRACTICAL_WORKFLOW.md) · [中文手册](USER_GUIDE.zh-CN.md)

PACE 根据实测调控活性和启动子接触证据，对目标基因的候选调控元件进行排序。分数表示候选背景内的相对支持，不是因果概率，也不是增强子贡献的基因表达比例。

## 一、主公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)}
{\displaystyle\sum_{e\in\mathcal E(G)} A_\star(e)\,\overline C(e,G)}
}
```

分子为当前元件的活性与综合接触之积，分母为该基因**完整计划候选集**的支持总和，默认包含启动子单元。只有全部支持可计算且分母大于零，才输出 `pace_score`。存在缺失时主分数为 NA，已知子集的份额另列为 `pace_score_conditional`。这里的完整仅指输入候选集，不代表已发现全部真实增强子。

| 符号 | 含义 | 对应字段 |
|---|---|---|
| E、G | 当前调控元件、目标基因 | element_id、gene_id |
| e、E(G) | 分母中的一个候选、完整计划候选集 | candidates.tsv |
| M、m | 固定检测组合、其中一种检测 | activity.panel |
| x_m(E) | 重复汇总后的实测信号 | resolved_activity.tsv |
| A_star(E) | 活性几何平均 | A_used |
| T(G)、t | 去重后的物理 TSS 集合、其中一个 TSS | promoters.tsv |
| pi(t given G) | 基因内总和为 1 的 TSS 权重 | promoter_weights.tsv |
| Ctilde(E,t) | 处理后的元件—TSS 接触 | resolved_contacts.tsv |
| Cbar(E,G) | 按所选 TSS 权重汇总的接触 | Cbar |

## 二、活性

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|}.
```

支持 ATAC、DNase、H3K27ac 单层，以及 ATAC+H3K27ac、DNase+H3K27ac 双层。所有元件使用同一组合；缺一层为 NA，实测零默认使活性为零。先在生物学重复内平均技术重复，再在供体内平均生物学重复，最后对供体等权平均。汇总前必须保证信号归一化可比；单个体使用 `individual`，多个体汇总使用 `population_mean`。

若某一检测在全部元件上的信号都乘以相同正数，全部活性会乘以一个共同系数，在同一基因的分子、分母中抵消。该性质不能消除信噪比、窗口、缺失模式或重复间尺度差异。

浅测序时可以显式添加分检测类型的伪计数：

```math
A_{\star,\epsilon}(E)=
\left[\prod_{m\in\mathcal M}(x_m(E)+\epsilon_m)\right]^{1/|\mathcal M|}.
```

`activity.pseudocounts` 默认为空，即全部 epsilon 为零。伪计数与对应归一化信号同单位，在重复汇总后加入；原始实测值保留，添加量写入 `activity_pseudocount`。缺失值始终为 NA。改变检测尺度时必须同时缩放伪计数，才能保持上述抵消性质。正伪计数可能增加背景支持，需比较关闭和启用时的结果，不宜套用统一数值。

`pace normalize-activity` 实现“窗口计数 × 10^6 / 过滤后文库片段数 / 窗口长度”。已经归一化的 bigWig 不应重复归一化；文库大小校正不等于批次校正。

## 三、接触

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}.
```

| 情形 | 使用的接触 |
|---|---|
| 同 bin 或指定近距离范围，有兼容先验 | P(d) |
| 同 bin、无先验，使用默认 prior_or_neighbor | cooler 接口记录的有效相邻接触最大值 |
| 近对角没有可用替代值 | NA |
| 非对角有限实测值，包括零 | 实测值加配置允许的伪计数 |
| 接触缺失，允许回退且有兼容先验 | P(d)，明确标注先验来源 |
| 接触缺失，不允许回退 | NA |
| prior_only | P(d) |
| shrinkage | r C_obs + (1-r) P(d) |

先执行近对角处理。同 bin 取决于分辨率和 bin 边界，不等于所有 ±5 kb 内的元件。`unresolved` 保留近对角缺失；`prior_or_unresolved` 只允许匹配先验；默认 `prior_or_neighbor` 还允许带来源记录的邻近最大值。邻近最大值是接触替代量，不是已验证的增强子—启动子环。

`contact.pseudocount: auto` 仅在 observed 模式且有兼容先验时，为非对角有限实测值添加 p(d)。`powerlaw` 强制要求先验，`none` 关闭添加。默认 kappa=1、d0=5000 bp，尚非经畜禽数据优化的参数。近对角替代值、prior_only 和 shrinkage 不再重复加伪计数。无效平衡 bin 不当作零，除非显式允许先验替代，否则保持不可用。

结果保留 `observed_value`、`prior_value`、`pseudocount_value`、来源、先验编号和处理原因。添加伪计数的接触标为 `regularized`，保留的观测系数不代表置信概率。

### 先验拟合和跨组织使用

`pace fit-prior --cooler ...` 从距离分箱均值拟合衰减，所有可测 bin 对的零计数也进入均值。默认拟合下限为矩阵分辨率，因此 10–25 kb 数据无需沿用 5 kb 下限。全零分箱不进入对数回归，但保留在报告中。保留染色体可用于检查接触衰减拟合，不能替代调控联系的功能验证。

先验与实测接触混用时，尺度、归一化、分辨率、平衡方式及窗口必须一致。其他组织先验默认拒绝；设置 `contact.allow_cross_context_prior: true` 后，可使用同物种、同组装、同目标层级的其他组织先验。保留来源组织、记录目标组织，并标注未经目标背景验证的迁移；`validated` 模式不接受这种迁移。仅修改组织名或尺度名不能完成校准。

完全没有 Hi-C 时，可显式选择 `mode: prior_only` 和 `prior_preset: abc_human`，仍须有实测活性。该选项采用 [ABC 官方配置](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml)的 gamma=1.024238616787792，并令 a=1、d_ref=d_min=5000，使用相对尺度。它不是通用畜禽参数，仅供 research 模式下的明确对照，不能与实测接触表混合。幅度 a 在仅先验评分中抵消。

## 四、多 TSS

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t).
```

相同物理 TSS 去重。同一实测 bin 内的不同 TSS 可复用相同元件、样本的接触查询，不增加重复数；按精确坐标计算的距离先验仍可能不同。默认要求全部正权重 TSS 的接触可用，缺一个则 Cbar 为 NA；零权重 TSS 不要求接触。

`pace prepare-promoter-weights` 可按 ATAC、DNase、H3K4me3 或 CAGE 启动子实测信号生成比例权重。这些权重是启动子使用的代理量；基因 TPM 本身不能确定 TSS 使用比例。

可选筛选始终使**同一基因的所有候选共用一套 TSS**：

- `minimum_weight`：删除原始权重低于阈值的 TSS，默认 0。
- `missing_policy: drop_missing`：某 TSS 对任一计划候选缺少可解析接触，就从该基因全部候选的接触汇总中删除。默认 `strict` 不删除。
- 剩余原始权重总和须达到 `minimum_retained_weight`，默认 0.9；不足时该基因不评分。

```math
\pi_{\mathrm{used}}(t\mid G)=
\frac{\pi(t\mid G)}{\sum_{u\in\mathcal T_{\mathrm{keep}}(G)}\pi(u\mid G)},
\qquad t\in\mathcal T_{\mathrm{keep}}(G).
```

候选元件和启动子单元仍保留在原分母目录。`promoter_weights.tsv` 记录原始及实际权重、删除原因、缺失候选数和剩余权重。`tss_contact_scope=selected_tss_set` 表示分数针对筛选后的启动子定义，不能解释为恢复了被删除启动子的调控；敏感性区间也仅针对该定义。不同筛选结果不能直接作为完整生物学差异进行 `pace compare`。对于很稀疏的接触图，统一删除可能无法保留足够权重，此时应使用有依据的先验或保留默认 NA 结果。

## 五、实验性跨基因分配

```math
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)},\qquad
\mathrm{PACE}_{\eta}(E,G)=
\frac{A_\star(E)\overline C(E,G)B(E,G)^\eta}
{\sum_{e\in\mathcal E(G)}A_\star(e)\overline C(e,G)B(e,G)^\eta}.
```

eta=0 时省略分配项及其数据要求，得到主公式。`allocation.eta: auto` 在没有合格独立功能标签时取零；显式非零 eta 属于实验性选择。eta 非零时，某元件任一候选目标基因的接触缺失都会使其 B 无法确定，不能静默删除该目标基因。

Cbar×B^eta 等于 Cbar^(1+eta)/(sum_H Cbar)^eta，既增强接触差异，也受到候选基因密度和注释影响。目前不能视为已证实的生物学竞争规律。应在独立功能数据上与 eta=0 及 `contact_power_2` 对照比较，见[校准说明](eta_calibration.md)。

## 六、分母缺失

令 S=A_star×Cbar，或用户显式选择的实验性扩展支持，则：

```math
\mathrm{PACE}_{\mathrm{conditional}}(E,G)=
\frac{S(E,G)}{\sum_{e\in\mathcal E^{\mathrm{score}}(G)}S(e,G)}.
```

`scoring.partial_policy: conditional` 可将条件分数也写入主列以兼容旧分析，但 `score_scope` 仍明确标注条件结果。默认 `withhold` 不这样做。例如已知支持为 2 和 1，条件分数为 2/3；若未测候选的实际支持是 97，完整分数仅为 0.02。

给定非负支持范围 L_i <= S_i <= U_i：

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}.
```

已解析支持固定，未知支持默认 [0,infinity)。可用 `inputs.support_bounds` 提供有依据的上下界及来源，不会据此填补主分数。支持为 (2,1,未知) 时，第一个元件的完整份额范围为 [0,2/3]；未知支持约束为 [2,4] 后为 [2/7,2/5]。这是总支持为正条件下的敏感性范围，**不是置信区间**。全零总支持输出 NA；未知项间的相关性可能使区间偏保守。

## 七、候选集与其他组学

网格只覆盖输入区域和必要启动子，不进行全基因组平铺。重复、重叠峰不会复制单元。过宽或噪声较多的输入区域仍会扩大分母，因此上游峰质量必须控制。`region_scores.tsv` 对区域包含的唯一单元求和，不重建分母；重叠区域不能再作为独立单元相加。该目录与 ABC 的峰顶窗口目录不同，不宜直接套用 ABC 的经验阈值。

RNA、其他组蛋白、CTCF 和甲基化默认作为注释，或用于独立分类器。相同基因表达权重同时乘入分子、分母会抵消；只在归一化后相乘则改变分数含义。

参考：[Fulco 等，2019](https://doi.org/10.1038/s41588-019-0538-0)、[Nasser 等，2021](https://doi.org/10.1038/s41586-021-03446-x)、[ABC 官方方法](https://abc-enhancer-gene-prediction.readthedocs.io/en/stable/usage/methods.html)。人类基准结果不能代替畜禽独立验证。
