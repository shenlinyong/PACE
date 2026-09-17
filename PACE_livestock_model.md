# PACE 家养动物调控预测模型

本文件是当前模型入口。完整定义统一维护于[数学公式](docs/FORMULA.md)，
生物学条件、三种模式、估计对象与限制见[模型说明](docs/model.md)。

$$
\boxed{
\operatorname{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
$$

默认没有适用功能校准时实际使用 eta=0；有适用训练/校准标签时估计
连续的 eta∈[0,1] 并冻结，测试集不参与估计。也可以预先指定固定参数。
分母仅为目标基因实际可评分候选的支持总和。零支持总和返回 NA。

实测、混合和基因组预测共用此公式。固定检测层的活性几何均值、
多TSS接触加权、候选集合与缺失规则均以当前规范为准。
[参数](docs/parameters.md) · [eta校准](docs/eta_calibration.md) ·
[命令行](docs/cli.md) · [独立手算](PACE_Review_and_Validation.md) ·
[实际验证状态](docs/validation.md)。

软件已提供三模式计算、适用资产检查、训练/校准和可复现测试；
仓库不附带已验证的真实家养动物权重，不声称软件测试证明生物学准确率。
