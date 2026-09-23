# PACE 总公式与详细解读

[英文公式](FORMULA.md) · [实际操作](PRACTICAL_WORKFLOW.md) · [完整中文手册](USER_GUIDE.zh-CN.md)

PACE 以**实测活性**为基础，计算一个候选调控元件对目标基因的相对支持。结果不是因果概率，也不是基因表达量中由该元件贡献的比例。

## 一、总公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

分子是当前元件的“活性 × 综合接触 × 可选分配项”；分母是该基因**完整计划候选集**的同类支持总和。默认只有候选支持全部可计算、且总支持大于零时，才输出主分数 `pace_score`。

如果分母有候选缺失，主分数为 NA；另外输出 `pace_score_conditional` 和 `pace_score_lo/hi`。这样不会把“只看到了部分候选”的高分当成完整背景下的高分。完整候选集仍由研究者预先定义，不等于全部真实增强子。

| 符号 | 含义 |
|---|---|
| E、G | 当前评分单元和目标基因 |
| e、E(G) | 逐个遍历的候选单元、该基因完整计划候选集 |
| A_star(E) | 固定检测组合的实测活性几何平均 |
| Cbar(E,G) | 对不同物理 TSS 的接触按固定权重汇总 |
| B(E,G) | 当前元件对该基因的接触占其全部候选基因接触的份额 |
| eta_used | 实际使用的可选分配指数；无合格独立功能证据时为零 |
| pi(t given G) | 预先固定的启动子权重，每个基因内加总为 1 |
| Ctilde(E,t) | 按明确策略校正或补充后的元件—TSS 接触 |

## 二、完全展开的总公式

```math
\mathrm{PACE}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(E,u)}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(e,u)}
\right]^{\eta_{\mathrm{used}}}
}.
```

每一个 Ctilde 都按下面相同的接触规则计算；每一个 x_obs 都来自真实实验。eta=0 时整个分配因子直接省略，不需要先计算 B。

## 三、活性如何计算

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}.
```

支持 ATAC、DNase、H3K27ac 单层，或 ATAC+H3K27ac、DNase+H3K27ac 双层。双层信号 4 和 9 的活性为 6。所选层缺失则为 NA；不能逐个元件临时删掉缺失层。

先汇总同一生物重复的技术重复，再汇总供体内生物重复，最后按明确的群体目标对供体等权。几何平均在这些测量汇总之后进行。多个个体要明确使用 population_mean。

原始窗口计数可用 `pace normalize-activity` 转为每百万合格片段、每碱基的信号。分母必须是完整合格文库的片段数，不能用候选峰内计数总和代替。已归一化 bigWig 不再重复做 CPM；文库量校正不能消除批次效应。

## 四、接触如何计算

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}.
```

P(d) 是距离背景，p(d) 是用于稀疏实测接触的伪计数。普通非对角实测位置采用：

```math
\widetilde C(E,t)=C_{\mathrm{obs}}(E,t)+p(d(E,t)),\qquad
\overline C(E,G)=\sum_t\pi(t\mid G)\widetilde C(E,t).
```

默认 `pseudocount: auto`：有匹配先验且使用 observed 接触模式时加伪计数，否则保留实测值并报告先验不可用。`none` 可关闭，`powerlaw` 强制要求匹配先验。默认 kappa=1、d0=5000 bp 是算法设置，尚不能称为畜禽最优参数。接触计数为零可能来自浅测序，不能据此判定没有生物学互作；原始零值仍保留在输出中。

| 数据情况 | 处理 |
|---|---|
| 同 bin 或指定近距离范围，且有匹配先验 | 使用 P(d)，记录 near_diagonal_prior，不再重复加伪计数 |
| 同 bin、无先验，但 cooler 提取时获得有效邻近接触 | 默认用有记录的邻近接触最大值，记录 near_diagonal_neighbor_max |
| 近对角没有可用校正依据 | NA；该基因主分数不再按不完整分母输出 |
| 普通位置有实测接触 | 按设定使用实测值加伪计数，或保留原值 |
| 普通位置缺接触且显式允许先验回退 | 使用匹配 P(d)，标明先验来源 |
| 普通位置缺接触且未允许回退 | NA |
| 显式 prior_only | 使用 P(d)，仍必须有实测活性 |
| 显式 shrinkage | 使用 r C_obs+(1-r)P(d)，不再叠加伪计数 |

近对角策略先于其他接触选择。它与 `allow_prior_fallback` 控制的普通位置回退不同。有匹配先验时，关闭普通回退并不会阻止近对角使用先验。同 bin 取决于实际分辨率和边界，不等于所有距离小于 5 kb 的联系。

加伪计数的证据类型为 regularized，区别于原始测量和凸组合收缩。`observed_value`、`prior_value`、`pseudocount_value`、模型身份、校正原因全部保留。这里 regularized 的 reliability=1 表示原观测的系数仍为 1，不表示后验可信度为 100%。

`pace fit-prior` 可以直接从 cool/mcool 拟合先验，将可测 bin 的未存储零计数计入距离分箱均值；无效平衡 bin 不参与。全零距离分箱不能进入对数回归，会在报告中列出。留出染色体只能评价接触背景拟合，不能替代增强子功能验证。先验必须与接触的尺度、分辨率、归一化、平衡方式和窗口一致。

没有 Hi-C 时，可**主动选择** `prior_preset: abc_human` 和 `mode: prior_only`。它使用人类 ABC 公布的 gamma=1.024238616787792，a=1、d_ref=d_min=5000，并明确标为 `abc_human_default` 和未验证迁移。其绝对尺度是相对量，不能混入实测 Hi-C；只能用于 research，不能当作已验证的家养动物参数。所有情况都保留实测活性要求，不恢复仅基因组预测。

## 五、多个启动子和可选分配项

```math
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
```

相同物理 TSS 先去重；不同 TSS 如果落在同一个实测 bin，可以复用同一元件、同一样本的可用接触查询，不增加重复数。保留各 TSS 的生物学身份和固定权重，另报告 `n_contact_bins`。距离先验使用真实距离，同 bin 内的不同 TSS 不一定有相同先验。

可通过 `pace prepare-promoter-weights` 用 ATAC、DNase、H3K4me3 或 CAGE 启动子信号生成比例权重。这是实验信号支持的代理权重，不自动等同于真实启动子使用率。所有计划 TSS 均需有测量；全零基因默认报错，只有明确选择时才退回等权。基因总 TPM 不能推断 TSS 使用比例。

B 的数学作用确实包含接触指数变化和候选基因背景校正：Cbar×B^eta=Cbar^(1+eta)/(sum_H Cbar)^eta。它会受注释版本与候选基因集合影响，因此默认无验证数据时关闭。基准新增 contact_power_2，对照单纯把接触平方的效果；B 的独立增益需要在冻结注释、训练/测试分离后证明，不能凭公式宣称。

## 六、缺失候选时的详细公式

```math
S_i=A_\star(i)\overline C(i,G)B(i,G)^{\eta_{\mathrm{used}}},\qquad
\mathrm{PACE}_{i,\mathrm{conditional}}=\frac{S_i}{\sum_{j\in\mathcal E^{\mathrm{score}}(G)}S_j}.
```

条件分数只在可用背景内归一化。默认主列为 NA，条件值在独立列中保留；显式设置 `scoring.partial_policy: conditional` 可以兼容旧的探索性用法，`score_scope` 仍会标记条件背景。

如果各候选支持满足 L_i≤S_i≤U_i，则：

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}.
```

已解析支持的上下界相等。未知支持默认 [0,无穷)，不会被假定成零。可以通过 support_bounds.tsv 提供有独立依据的**最终支持值**上下界及 bound_source；这不会填补实测活性或使主分数变为完整。上界 NA 表示无界。

例如已知支持 2、1，另一个未知：第一个条件分数是 2/3，但完整背景下只能给 [0,2/3]。若外部依据将未知支持约束在 [2,4]，区间变为 [2/7,2/5]。没有依据时不能自动制造窄区间。

这些是以支持范围假设为条件的**敏感性范围，不是统计置信区间**。候选间有依赖时可以偏保守。区间以总支持大于零为条件；全部可能支持均为零时返回 NA，唯一可能有支持的候选区间为 [1,1]，也不代表获得了生物学验证。

## 七、峰区域与其他组学

一个峰可能跨多个唯一网格单元。`region_scores.tsv` 对同一来源、区域和基因的单元只加总一次，不创建第二套分母；重叠区域之间不能再作为独立元件求和。区域区间由单元区间保守加总，上界最多为 1。这种网格方法不等同于 ABC 以峰顶为中心的候选窗口。

RNA、其他组蛋白、CTCF、WGBS/RRBS 仍有接口，默认作注释或独立验证的分类器特征；不任意乘入抑制系数或表达权重。软件改进能改善可追溯性和计算边界，是否提高生物学准确率仍需按[验证方案](ROBUSTNESS_VALIDATION.md)检验。
