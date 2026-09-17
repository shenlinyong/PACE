# PACE

PACE 是用于研究家养动物增强子—基因调控联系的 Python 软件。它把活动信号、启动子接触和可选的跨基因分配组合为**基因内部的相对增强支持份额**，保留每条候选边的来源、缺失原因和实际分母。

作者与维护者：**申林用（Linyong Shen，shenlinyong），西北农林科技大学**。

## 当前总公式

$$
\boxed{
\operatorname{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
$$

分母只有目标基因实际可评分候选的支持总和。默认无适用功能校准时 eta=0；
有合格训练/校准数据时估计 [0,1] 内的连续值。测试集不参与。
完整定义见[数学公式](docs/FORMULA.md)。

## 安装与最小示例

需要 Python 3.11 或更新版本。以下命令从仓库安装；暂不表示已经发布到 PyPI。

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
python -m venv .venv
source .venv/bin/activate
python -m pip install .
PACE --version
PACE --mode measured --config examples/measured/config.yaml --out results/measured
PACE --mode hybrid --config examples/hybrid/config.yaml --out results/hybrid
PACE --mode genome --config examples/genome_only/config.yaml --out results/genome
```

三个示例均离线运行，使用明确标记的合成数据和固定测试权重。基础安装不需要 GPU。
读取 bigWig、cool/mcool、BCF 使用 `pip install '.[io]'`；CNN 训练和推理使用
`pip install '.[sequence]'`。安装依赖可能需要网络，运行示例不下载数据或模型。

```bash
pace-livestock validate --config examples/measured/config.yaml
pace-livestock run --config examples/measured/config.yaml --out results/measured
pace-livestock capabilities --config examples/genome_only/config.yaml
```

也可以用 `./install.sh --prefix "$HOME/.local" --python python3` 安装独立环境。
`PACE measured`、`PACE hybrid`、`PACE genome` 与 `--mode` 写法等价。
YAML 可省略，直接指定 `--catalog-dir`、`--activity`、`--contacts`、`--reference`、
`--sequence-model` 等参数；完整真实数据示例见[命令行手册](docs/cli.md)。
`PACE run --help` 显示全部参数，原 `pace-livestock` 命令继续可用。

YAML 内路径以配置文件所在目录为基准，命令行路径以当前目录为基准；
命令行参数覆盖配置值，已有输出目录不会被静默覆盖。

## 三种模式

| 模式 | 输入与用途 |
|---|---|
| measured | 规范单元、候选边、真实样本活动表、接触表或有来源的接触先验 |
| hybrid | 合格实测与适用的定量序列预测；有匹配校准器才融合，否则按明确规则选择单来源 |
| genome_only | atlas/启动子候选、参考 FASTA、适用定量权重和接触先验；个体化另需 VCF/BCF、可调用区间、倍性及可靠相位 |

主流程统一为：逐检测层得到 bulk 边际信号 → 活性几何均值 → 多 TSS 接触 → 可选分配 → 基因内归一化。
默认 `--eta auto`：没有适用功能标签时实际 `eta=0`；通过 `--eta-labels` 提供标签时，
从训练/校准数据自动估计连续 `0≤eta≤1`，测试集不参与估计。`--eta-model` 可冻结复用，
也可用 `--eta 0.35` 等数值预先指定。详见[公式与校准说明](docs/eta_calibration.md)。
H3K27ac-only 等单层模式可以显式配置，双层运行中缺少一层不会自动变成单层。
RNA、H3K4me1、H3K4me3、H3K27me3、CTCF 和甲基化默认是注释；有功能标签时可以进入独立分类器。

## 读懂结果

- `scores.tsv.gz`：全部候选边、A/C/B、原尺度及对数支持、分母、PACE、缺失与来源。
- `gene_summary.tsv`：每个基因的候选数、可评分数、条目覆盖与实际归一化集合 ID。
- `resolved_activity.tsv` / `resolved_contacts.tsv`：实测、预测、融合或先验的解析结果。
- `multiomics_features.tsv.gz`：附加组学注释及角色。
- `qc_report.json` / `report.md`：能力、覆盖、不可支持变异和解释边界。
- `run_manifest.json` / `resolved_config.yaml`：输入与模型哈希、环境、参数和随机种子。
- `eta_calibration.json`：实际 eta、拟合/回退原因、标签来源、适用范围及冻结复用记录。

技术缺失是 NA，真实零是 0。`partial` 分数只针对可测子集；分数和为 1 不表示数据完整。
PACE 不是因果概率，也不直接预测表达变化。增强子自身支持不变时，其他候选变化仍可改变它的份额。
跨样本比较使用 `compare` 从 support 重算共同分母，分别给出完整与条件 Delta。

## 模型资产与适用边界

本仓库提供真实训练和推理代码、合成示例及软件测试，**不附带已验证的真实家养动物权重**。
不能把 demo 成功解释为已证明跨品种、跨组织或个体变异预测准确。

基因组流程支持 SNV 和保持固定目标对应关系的短 indel。缺失 GT、无法连接的相位区块、不可调用区域、
目标被 indel 改变或已报告但无法解析的 SV，会返回明确不可用状态。软件不承担比对、变异发现、SV calling、
全基因组新增强子发现或复杂 CNV 剂量重建。单变异 REF/ALT 情景与完整个体效应分开输出。

当前分支只有一套评分实现，所有命令别名和 `scripts/pace.py` 都调用相同的软件包。
已淘汰的工作流和公式仅通过 Git 历史追溯，参见[迁移说明](docs/migration.md)。

详细说明：[新软件手册](docs/software.md) · [输入准备](docs/input_preparation.md) · [训练与校准](docs/training.md) ·
[比较与基准](docs/comparison.md) · [测试记录](docs/validation.md) · [数学合同](docs/model.md)。

引用软件时请记录 Git commit；作者和引用信息见 [CITATION.cff](CITATION.cff)。当前没有软件论文 DOI。
