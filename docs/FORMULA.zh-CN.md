# PACE 公式与原理：中文说明

本文解释当前软件实际计算的模型，并给出通用总公式和三种数据条件下**分别展开的详细总公式**。数学表达式与[英文数学定义](FORMULA.md)完全一致。安装、配置和命令见[完整中文使用说明](USER_GUIDE.zh-CN.md)。

## 1. 这个模型到底在算什么

对一个给定基因，软件先计算每个候选元件的调控支持，再求它占这个基因全部可评分候选支持的比例。

例如，同一基因有两个候选元件，支持分别为 12 和 2，则分数分别为 12/14≈0.857 和 2/14≈0.143。前者在这个背景中的支持更大；这不表示它有 85.7% 的因果概率，也不表示它贡献了 85.7% 的基因表达。

评分对象是固定坐标的非重叠单元。一个单元可以对应增强子、启动子或二者共有的区域。默认纳入合格启动子单元，因此分母应准确称为“计划候选调控单元”，不能宣称覆盖了全部生物学增强子。

## 2. 简写总公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

分子的未归一化支持记为：

```math
S(E,G)=A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}.
```

依次理解三个因子：

1. **A_star：元件有多活跃。** 从本次固定的检测层组合估计活性。
2. **Cbar：它与目标基因的启动子有多少接触支持。** 可以综合同一个基因多个物理 TSS。
3. **B 的指数项：可选的目标分配修正。** 描述该元件的接触有多大份额偏向这个基因。

最后除以同一基因可评分候选的支持总和。默认自动设置在证据不足时使用数值 0，因此可选分配项整体跳过；有合格功能数据通过分组验证后才使用 [0,1] 内的学习结果。这里的核心仍是活性、接触和明确的归一化背景。

## 3. 所有符号分别代表什么

| 符号 | 中文含义 | 对应输入或输出 |
|---|---|---|
| E | 当前评分元件 | element_id |
| G | 当前目标基因 | gene_id |
| e | 分母中逐个遍历的其他候选元件，也包括 E | 该基因可评分集合 |
| H | 元件 E 的另一个候选目标基因 | candidates.tsv |
| M、m | 本次固定活性组合、其中一项检测 | activity.panel |
| x_obs,m | 检测层 m 的合格实测信号 | observed_activity 解析结果 |
| x_hat,m | 序列模型预测的定量信号 | predictions 或本地序列模型 |
| x_star,m | 最终选定或校准融合后的信号 | resolved_activity.resolved_value |
| A_star | 对所有必需检测层的最终信号取几何平均 | A_used |
| T(G)、t | 基因 G 的不同物理 TSS 集合、其中一个 TSS | promoters.tsv |
| pi(t given G) | 这个 TSS 在该基因接触汇总中的固定权重 | pi；同一基因内和为 1 |
| C_obs | 实测接触 | observed_contacts |
| C_prior | 适用接触先验 | contact prior 资产 |
| r | 接触组合中实测来源的权重 | reliability；范围 [0,1] |
| Cbar | 经来源处理并汇总多 TSS 的基因接触 | Cbar |
| G(E) | 预先固定的元件 E 候选基因集合 | candidates.tsv |
| B | E 对 G 的接触占 E 对所有候选基因接触的比例 | B；实际使用时输出 |
| eta_used | 本次真正使用的分配指数 | eta_calibration.json |
| E_score(G) | 预定候选中本次能够计算支持的集合 | scoreable 和实际分母身份 |
| S | 归一化前的支持 | support / log_support |

**A-hat 与 A-star 的区别：**帽子表示“预测”，星号表示“最后用于评分的值”。仅当最终选择序列预测时二者一致。软件先在每个检测层解析来源，再构造活性；不能把这两个符号理解为同一量的不同写法。

## 4. 活性是怎样得到的

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\star,m}(E)\right]^{1/|\mathcal M|}.
```

若组合是 ATAC+H3K27ac，就是两者乘积开平方；若只选择 H3K27ac，就直接使用 H3K27ac 信号。支持的组合为 ATAC、DNase、H3K27ac 单层，或 ATAC+H3K27ac、DNase+H3K27ac 双层。

整个运行保持同一个组合。双层分析中某个元件缺少 ATAC，不能只给这个元件改成 H3K27ac 单层。实测零保持零，未测或不合格保持 NA。技术重复、生物重复和供体分别聚合，避免技术重复较多的动物获得更大权重；individual 目标不能混合不同供体。

混合模式有合适校准器时，在每一层采用以下公式：

```math
x_{\star,m}(E)=s_m\left\{
\exp\left[
\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(E)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(E)}{s_m}\right)
\right]-1\right\},\qquad s_m>0,\quad 0\le\lambda_m\le1.
```

s_m 是校准时冻结的正尺度；lambda_m 是该检测层、该质量分组的实测权重。它们必须与当前组织、单位、归一化和窗口匹配。这里的 log1p 变换只用于**信号融合**，不是向主活性或分母添加任意常数。

没有可识别的适用校准器时，优先使用合格实测；实测不可用且有适用预测时使用预测。已校准的内部权重需要两种来源，若其中一种缺失则按上述单来源规则回退；明确的端点权重只要求被选中的来源。两种来源都不可用时保留 NA。软件不会自行指定 50:50 融合。

## 5. 接触、多启动子和分配项

首先处理每个元件—TSS 的接触来源：

```math
\widetilde C(E,t)=r(E,t)C_{\mathrm{obs}}(E,t)
+[1-r(E,t)]C_{\mathrm{prior}}(E,t),
```

r=1 表示只用实测，r=0 表示只用先验，介于两者之间表示按已声明的可靠性收缩。实测与先验必须具有兼容的分辨率、量纲、归一化和背景，不能把两个不相容的数字直接平均。软件只要求实际有正权重的来源，不计算 0×NA。

再汇总该基因的各个物理 TSS：

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),
\qquad \sum_{t\in\mathcal T(G)}\pi(t\mid G)=1,
```

pi 是 TSS 权重，与实测可靠性 r 不同。多个转录本共用一个 TSS 时只计一次；没有可信起始使用比例时可以显式选择等权。只知道基因 TPM 不能反推出各 TSS 使用率。某个正权重 TSS 的接触缺失时，不把它删除后重新分配权重。

可选分配项为：

```math
B(E,G)=\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\displaystyle\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)\widetilde C(E,u)}.
```

这个分母跨的是**同一元件的候选基因**；PACE 总公式的分母跨的是**同一基因的候选元件**。两种归一化回答不同问题，不能混淆。B 是分配指标，不代表物理资源守恒。有正分配指数时，必须保留预定候选基因集合，不能通过删除缺接触的基因抬高 B。

## 6. 完全展开的通用总公式

以下把活性几何平均、多 TSS 接触和跨基因分配全部代入，不再只写 A_star、Cbar 和 B 的缩写。

```math
\mathrm{PACE}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\star,m}(E)\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]\right\}
\left\{
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)
[r(E,u)C_{\mathrm{obs}}(E,u)+(1-r(E,u))C_{\mathrm{prior}}(E,u)]}
\right\}^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\star,m}(e)\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]\right\}
\left\{
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)
[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)
[r(e,u)C_{\mathrm{obs}}(e,u)+(1-r(e,u))C_{\mathrm{prior}}(e,u)]}
\right\}^{\eta_{\mathrm{used}}}
}.
```

公式虽然长，计算过程没有增加：逐层得到 x_star → 计算活性 → 处理并汇总接触 → 按需要计算分配 → 对基因归一化。x_star 的来源由下列三种情况确定。

## 7. 情况一：实测数据充分，measured

以下表示活性和接触均由合格实测支持的情况，即 x_star=x_obs、r=1：

```math
\mathrm{PACE}_{\mathrm{measured}}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(E,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(E,t)}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)C_{\mathrm{obs}}(E,u)}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(e,t)\right]
\left[
\frac{\sum_{t\in\mathcal T(G)}\pi(t\mid G)C_{\mathrm{obs}}(e,t)}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)C_{\mathrm{obs}}(e,u)}
\right]^{\eta_{\mathrm{used}}}
}.
```

这里没有要求序列模型。只有 H3K27ac 时可在整个运行中明确设为单层；拥有 ATAC 或 DNase 与 H3K27ac 时可使用双层几何平均。

measured 还允许“实测活性+明确先验接触”或接触收缩。此时仅把上式的实测接触部分替换为通用总公式中的实测/先验组合，活性仍由实测给出。近对角线接触按显式先验或缺失策略处理，不自动视为可信互作。

## 8. 情况二：实测不完整或质量有限，hybrid

当每个必需活性层都成功融合，详细总公式如下：

```math
\mathrm{PACE}_{\mathrm{hybrid}}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}s_m\left\{
\exp\left[\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(E)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(E)}{s_m}\right)\right]-1\right\}\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]\right\}
\left\{\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(E,t)C_{\mathrm{obs}}(E,t)+(1-r(E,t))C_{\mathrm{prior}}(E,t)]}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)[r(E,u)C_{\mathrm{obs}}(E,u)+(1-r(E,u))C_{\mathrm{prior}}(E,u)]}
\right\}^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}s_m\left\{
\exp\left[\lambda_m\log\left(1+\frac{x_{\mathrm{obs},m}(e)}{s_m}\right)
+(1-\lambda_m)\log\left(1+\frac{\widehat x_m(e)}{s_m}\right)\right]-1\right\}\right]^{1/|\mathcal M|}
\left\{\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]\right\}
\left\{\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)[r(e,t)C_{\mathrm{obs}}(e,t)+(1-r(e,t))C_{\mathrm{prior}}(e,t)]}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)[r(e,u)C_{\mathrm{obs}}(e,u)+(1-r(e,u))C_{\mathrm{prior}}(e,u)]}
\right\}^{\eta_{\mathrm{used}}}
}.
```

与 measured 相比，主要变化在活性来源：x_star 被替换为有校准依据的实测—预测融合表达式。接触可以来自实测、先验或显式可靠性组合。

若某一层按规则回退到单来源，就将对应乘积中的融合表达式替换为该层合格 x_obs 或 x_hat。不能为了填满结果而融合错组织、污染或尺度不相容的数据。这里的权重来自校准协议，不是软件自动创造的通用数据质量分数。

## 9. 情况三：新个体只有基因组数据，genome_only

当前实现使用定量序列模型预测各检测层，并采用适用距离先验计算接触：

```math
\widehat x_m(E)=\frac1{K_E}\sum_{k=1}^{K_E}f_{\theta,m}(\mathrm{seq}_{E,k}),
\qquad C_{\mathrm{prior}}(E,t)=a
\left[\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma}.
```

| 新符号 | 含义 |
|---|---|
| K_E | 元件所在位置的实际染色体拷贝数，当前个体路径支持 1 或 2 |
| seq_E,k | 第 k 个拷贝的合格固定目标序列窗口 |
| f_theta,m | 使用已训练参数 theta 预测检测层 m 定量信号的模型 |
| d_Et | 元件锚点到 TSS 的距离 |
| a | 距离先验的接触尺度参数 |
| gamma | 接触随距离下降的指数，要求为正 |
| d_min | 明确声明的近距离平台下限 |
| d_ref | 距离先验的参考距离 |

先对各拷贝**同一检测层的信号取平均**，然后计算活性；不能改成先算每个拷贝的支持再相加。序列模型输出必须是定量信号，峰概率或分类 logit 不能直接替代。距离先验参数需要匹配的资产，不内置所谓适合所有畜禽的人类常数。

详细总公式为：

```math
\mathrm{PACE}_{\mathrm{genome}}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}\left\{\frac1{K_E}\sum_{k=1}^{K_E}f_{\theta,m}(\mathrm{seq}_{E,k})\right\}\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}\right]
\left[\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{Et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
{\sum_{H\in\mathcal G(E)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)a\left(\frac{\max(d_{Eu},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
\right]^{\eta_{\mathrm{used}}}
}{
\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}\left\{\frac1{K_e}\sum_{k=1}^{K_e}f_{\theta,m}(\mathrm{seq}_{e,k})\right\}\right]^{1/|\mathcal M|}
\left[\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}\right]
\left[\frac{
\sum_{t\in\mathcal T(G)}\pi(t\mid G)a\left(\frac{\max(d_{et},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
{\sum_{H\in\mathcal G(e)}\sum_{u\in\mathcal T(H)}\pi(u\mid H)a\left(\frac{\max(d_{eu},d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}}
\right]^{\eta_{\mathrm{used}}}
}.
```

它表示**指定组织背景下的预测调控潜能**，并不测量这个个体当前的染色质状态。当前软件实现的是距离接触先验，不是序列预测三维互作的模型。

新个体缺表观数据时，可以在适用范围内使用已有同物种模型。整个物种没有功能训练数据时，则不能假定存在可用模型；跨物种迁移需要单独验证。只有 DNA 也无法唯一恢复营养、感染、激素或细胞比例等状态。

个体输入需明确参考、实际基因型、callability、倍性和可信相位。改变输出目标或 E–TSS 几何关系的变异、超出能力的 SV/BND 不会被假装成完整个体效应。导入外部预测仍要通过输入绑定和结构有效性检查。

## 10. RNA、甲基化和其他表观信号为什么单独处理

RNA-seq、H3K4me1、H3K4me3、H3K27me3、H3K9me3、CTCF 和 WGBS/RRBS 具有[明确软件接口](MULTIOMICS.md)。它们默认成为可追溯注释，有适用功能标签时可作为独立分类器的特征。

这些信号的作用不同：H3K4me3 多用于启动子状态，CTCF 与染色质结构有关，甲基化需要结合区域、覆盖及实验类型解释。给它们统一设置正负乘子会引入尚未验证的生物学假设，因此没有强行加入主支持公式。

同一基因的表达权重若同时乘进分子和分母，会完全抵消；若只在归一化之后相乘，则改变分数含义和跨基因尺度。当前软件保留清楚的主 PACE 份额，把 RNA 用于独立注释或预先声明的 ML 消融分析。独立 ML 的分数和校准概率不与 PACE 简单平均。

## 11. 零、缺失、分母和比较

| 情况 | 软件处理 | 正确解释 |
|---|---|---|
| 必需证据缺失或无效 | NA 支持及原因 | 不能算，不等于没有作用 |
| 必需证据明确且支持为真实零 | 保留零，并计入可评分集合 | 已知零与技术缺失不同 |
| 正分配指数需要的候选基因接触缺失 | 不可评分 | 不缩小 B 的候选集合 |
| 基因支持总和为零 | PACE 为 NA | 不加常数制造相对分数 |
| 部分计划候选不可评分 | normalization_status=partial | 仅对可评分子集归一化 |
| 所有计划候选可评分 | normalization_status=complete | 不等于找到了全部真实增强子 |
| 原尺度乘积极小或极大 | 保留 log_support 并在对数空间归一化 | 避免数值下溢/溢出破坏份额 |

前面的两个候选例子中，E1 活性为 sqrt(4×9)=6、接触为 2，支持为 12；E2 活性为 sqrt(1×4)=2、接触为 1，支持为 2。默认分配项关闭时，分数为 0.857 和 0.143。

如果 E2 后来因为技术缺失不能评分，E1 的条件分数会变成 1，但实验支持仍是原来的 12。跨样本必须检查候选覆盖并使用[共同背景比较](comparison.md)，同时报告活性、support、基因总支持及 PACE 份额变化。相对份额下降不能单独证明增强子活性或基因表达下降。

## 12. 对家养动物的适配体现在哪里

模型通过明确的物种/组装/组织背景、固定检测组合、实测与预测来源区分、多 TSS、供体重复、接触测量合同、个体倍性与基因型检查处理异质数据。对于表观检测少、参考版本多、启动子注释不齐、鸡性染色体倍性不同等情况，软件有对应的输入要求和失败解释。

这些设计使分析更可重复和可审查。它们不自动证明在所有家养动物、品种或组织上优于现有方法。软件回归测试检查计算行为；真实生物学性能、跨品种泛化和个体效应准确性仍需分别验证。详见[完整使用说明](USER_GUIDE.zh-CN.md)和[能力边界](limitations.md)。
