# PACE 实测活性模型：公式与解读

PACE 根据实测调控活性和明确来源的接触证据，对增强子等候选调控单元与基因的联系进行排序。
每次分析都需要 ATAC、DNase 或 H3K27ac 实测信号。主分数回答：在该基因的可评分候选背景中，这个元件占多少相对支持。

## 1. 总公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

先计算当前元件的“活性 × 综合接触 × 可选分配项”，再除以同一基因所有可评分候选的同类支持总和。
例如两条候选支持为 12 和 2，分数分别为 12/14 和 2/14。它们表示相对支持份额，不能解释成因果概率或表达贡献比例。

| 符号 | 含义 |
|---|---|
| E、G | 当前候选调控单元、目标基因 |
| e | 逐个遍历该基因的候选调控单元，包含 E |
| M、m | 本次固定的活性检测组合、其中一个检测层 |
| x_obs,m | 通过质控、统一量纲并按重复规则汇总后的实测信号 |
| A_star | 最终用于评分的实测活性，对应 A_used |
| T(G)、t | 基因 G 去重后的物理 TSS 集合、其中一个 TSS |
| pi(t given G) | 固定的 TSS 权重，同一基因内合计为 1 |
| Cbar | 按 TSS 权重汇总后的接触支持 |
| G(E)、H | 元件 E 的预定候选靶基因集合、其中一个基因 |
| B | E 对 G 的接触占 E 对所有候选靶基因接触的比例 |
| eta_used | 实际使用的分配指数，范围 [0,1] |
| E_score(G) | 本次能够计算支持的候选单元集合 |

## 2. 活性完全来自实测数据

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}.
```

可选组合是 ATAC、DNase、H3K27ac 单层，或 ATAC+H3K27ac、DNase+H3K27ac 双层。
例如 ATAC=4、H3K27ac=9，活性为 6。只有 H3K27ac 时，可以为整个运行明确选择 H3K27ac 单层。
已经选择双层后，某个元件缺少一层，其活性保留 NA，不能只在该元件上临时换成单层。

技术重复先在生物学重复内平均，再在供体内平均；群体汇总时供体等权。单个体分析不能合并多只动物的测量。
这些计算要求上游已完成适当归一化和质控，不自动消除批次效应。真实零保留为零，未测和不合格信号保留缺失。

## 3. 接触及多启动子汇总

默认采用实测接触。缺少合适 Hi-C 时，可显式提供适用的距离先验；也可在明确可靠性来源后进行接触收缩：

```math
\widetilde C(E,t)=rC_{\mathrm{obs}}(E,t)+(1-r)C_{\mathrm{prior}}(E,t),
\qquad
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t).
```

r=1 只要求实测接触，r=0 只要求接触先验，0<r<1 要求两种来源都有效。当前 r 是本次运行声明的权重，需要记录其估计或校准依据。
pi 是同一基因各 TSS 的权重，和 r 的用途不同。共用同一个物理 TSS 的转录本只计一次；没有可靠使用比例时可明确选择等权。

距离先验为：

```math
C_{\mathrm{prior}}(E,t)=a\left[\frac{\max(d(E,t),d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma}.
```

a、gamma、d_min 和 d_ref 必须来自适用的先验资产，不使用一套未经验证的通用畜禽常数。
采用先验时，活性仍必须实测，结果须注明接触来自距离先验。先验接触不等同于观察到了染色质互作。
实测接触、先验和导入结果的分辨率、量纲、归一化和窗口定义需要一致。

正权重 TSS 的接触缺失时，Cbar 保留 NA；不能删掉该 TSS 后重新分配权重。同一 Hi-C bin 内和近对角线的接触按显式先验或缺失策略处理。

## 4. 可选分配项 B

```math
B(E,G)=\frac{\overline C(E,G)}{\displaystyle\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
```

它比较同一元件的不同候选基因。主 PACE 分母则比较同一基因的不同候选元件。
默认自动设置在没有合格功能标签时使用 eta=0，直接跳过 B；有合格标签且通过独立分组验证后，才启用 [0,1] 内的学习结果。
使用非零指数时，必须保留预先确定的候选基因集合，不能通过删除缺接触的基因抬高 B。

## 5. 活性和接触均为实测时的完整展开式

```math
\mathrm{PACE}(E,G)=
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

这里已展开活性几何平均、多 TSS 汇总和跨基因接触分配。若显式选择接触先验或收缩，把式中每个 C_obs(E,t) 换为：

```math
rC_{\mathrm{obs}}(E,t)+(1-r)a\left[\frac{\max(d(E,t),d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma}.
```

分子、分母以及 B 内所有元件—TSS 接触都进行同样替换。x_obs 始终保持实测信号。

默认 eta=0 时，实际计算简化为：

```math
\mathrm{PACE}_{\eta=0}(E,G)=
\frac{
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(E)\right]^{1/|\mathcal M|}
\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t)}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
\left[\prod_{m\in\mathcal M}x_{\mathrm{obs},m}(e)\right]^{1/|\mathcal M|}
\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(e,t)}.
```

## 6. 其他组学及结果边界

RNA-seq、H3K4me1/3、H3K27me3、H3K9me3、CTCF、WGBS 和 RRBS 接口均保留，默认用于注释。
有适用功能标签时可以用于单独训练、验证的分类器，输出另列。没有固定表达乘子或通用甲基化惩罚系数；不把分类器分数与 PACE 相乘或平均。

完整候选表始终保留；只有部分候选可评分时，报告 partial 和实际分母。总支持为零返回 NA。
某个基因只剩一个正支持候选时，其分数可以为 1，但这不能证明该候选已被功能验证。
比较样本需要共同候选背景和可比较的检测条件，并重算共同分母。分数份额的变化不能单独证明基因表达或绝对活性的变化。

家养动物适配主要体现在允许固定单层活性、显式接触来源、多 TSS 去重、供体与重复管理、参考版本检查和缺失报告。
这些设计便于处理不齐全的实验数据，真实物种和组织中的预测表现仍需独立评估。

[完整使用说明](USER_GUIDE.zh-CN.md) · [多组学接口](MULTIOMICS.md) · [参数](parameters.md) · [数值例子](WORKED_EXAMPLES.md)
