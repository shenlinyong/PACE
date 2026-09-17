> **Interface scope:** This retained guide documents the legacy region-based scripts/workflow.
> For the installable canonical-grid package and its September 2026 defaults, see [software.md](software.md).

# PACE（Prediction of Activity-based regulatory Connections for Enhancers）

**结合增强子活性、启动子接触和基因注释，预测增强子–基因调控联系。**

[English](../README.md) · [安装依赖](INSTALLATION.md) · [完整参数](PARAMETERS.md) · [逐步教程](TUTORIAL.md) · [输入格式](INPUTS.md) · [结果字段](IO_FORMATS.md)

PACE 是可用于多物种的增强子–基因联系预测框架。它整合可用的活性测量，按可靠性调整接触信息，对同一基因的不同转录起始位点（TSS）进行去重和加权汇总，并分别报告预测分数与输入证据状态。

这些处理使 PACE **更适配家养动物中常见的数据条件**：多组学覆盖不齐、同组织同条件的 Hi-C 有限，以及启动子和候选目录注释不完整。其他物种在相同数据条件下也可使用这些功能。猪、牛和鸡是应用场景，模型本身不限定物种；运行时提供对应的参考基因组、注释和信号数据即可。

## 与原始 ABC 的直接比较

比较对象为 [Fulco 等（2019）的 ABC 模型](https://doi.org/10.1038/s41588-019-0538-0)及其 [NG2019 实现](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019)。ABC 为活性和接触能够可靠估计的场景提供直接的评分基线；PACE 在这一基础上，明确处理测量覆盖、可靠性和注释完整性的不均一。

| 原始 ABC 的处理 | PACE 的处理 | 修改原因及对家养动物数据的意义 |
| --- | --- | --- |
| 将选定的接触估计直接用于活性 × 接触；原研究已考虑平均 Hi-C 和距离替代 | 使用可靠性调整后的 $\overline C(E,G)$ | 匹配的 Hi-C 缺失或质量较低时，减少观测对评分的影响，并保留距离先验及接触来源 |
| 原始双测定方案对可及性与 H3K27ac 取几何平均 | $A(E)$ 按有效测量整合，明确区分 `NA` 与实测 `0` | 不同物种、组织和品种的测定覆盖不齐时，可保留已有测量；未测不再被当作无活性 |
| 对每个基因独立归一化 | 先加入 $B(E,G)$，在一个增强子的多个候选基因之间分配支持 | 降低广泛共享增强子的原始支持，旨在改善基因密集区域中靶基因判定的特异性 |
| 相对于选定的基因启动子计算接触 | 对不同 TSS 去重，以 $\pi(G,t)$ 加权平均 | 保留可变启动子信息，避免因重复转录本记录而放大支持 |
| 相对分数概括活性–接触支持 | 分数与 $Q$、证据状态及原因分别报告 | 注释、活性或接触证据不足时仍可保留候选，并标记为 `provisional` |
| 对候选活性–接触乘积归一化 | 定义 $\mathcal E^{\mathrm{obs}}(G)$，保留未评分行，并允许独立输入 $U(G)$ | 无法计算的候选仍是缺失值；候选目录外的支持是否已知也能被明确记录 |

“更适配”指模型和输出能够处理这些数据条件。预测准确率是否提高，需要在匹配数据预算、候选集合和验证规则的条件下独立比较。

原始 ABC 已包含接触归一化、距离伪计数、缺失值输出及失败基因记录。这里比较的是 PACE 新增的可靠性混合、有限支持分母、未评分计数和证据状态规则。[详细对照](ABC_COMPARISON.md#numerical-defaults-that-require-care-in-a-comparison)还列出两版距离指数、接触尺度和伪计数含义的差别，避免只因参数名称相近就沿用设置。

### 两个评分公式

$$
\mathrm{ABC}(E,G)=\frac{A_{\mathrm{ABC}}(E)C_{\mathrm{ABC}}(E,G)}{\sum_{e\in\mathcal E(G)}A_{\mathrm{ABC}}(e)C_{\mathrm{ABC}}(e,G)}.
$$

$$
\mathit{PACE}(E,G)=\frac{A(E)\overline C(E,G)B(E,G)^\eta}{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)\overline C(e,G)B(e,G)^\eta+U(G)}.
$$

这里的 $\overline C(E,G)$ 就是[正式公式](FORMULA.md)中的 $C(E,G)$、结果表中的 `contact_gene`：接触可靠性调整和 TSS 汇总均已完成。横线只用于与 ABC 的接触量区分，不代表另一层计算。

### 每项修改具体起什么作用

**接触可靠性。** 对增强子–TSS 对，计算 $C_{\mathrm{adj}}=P(d)[1-\lambda+\lambda H/D]$。$P$ 为距离先验，$H$ 为观测接触，$D$ 为兼容的期望接触，$\lambda$ 为独立评估的可靠性。当 $P=1$、$H/D=4$ 时，$\lambda=0、0.25、1$ 分别得到接触值 1、1.75、4。缺少合格观测时回退到先验；如果是可靠的实测零，接触可以降至零。

**有效活性。** 对预先缩放的信号做 `log1p` 加权平均，再用 `expm1` 转回。等权、质量为 1、尺度为 1 时，两个测定分别为 `(8, NA)` 得到活性 8；`(8, 0)` 得到 2；`(0, 0)` 得到 0；`(NA, NA)` 保持缺失。因此，缺失和实测零进入计算的方式不同。平移几何平均也改变了低信号行为，其假阳性影响需要验证。

**增强子侧分配。** $B(E,G)=\overline C(E,G)/\sum_g\overline C(E,g)$。若一个增强子对两个基因的接触为 3 和 1，则分配比例为 0.75 和 0.25。在活性为 1、$\eta=1$ 时，原始支持由 3 和 1 变成 2.25 和 0.25，再进行基因侧归一化。它并不保证所有最终分数都降低：如果每个基因都只有一个可评分候选且 $U=0$，两个分数仍可同时为 1。

**多 TSS。** 两个不同 TSS 的接触为 2 和 6，均匀加权得到 4；重复其中一个转录本记录不会改变结果。若独立启动子使用证据给出权重 0.75 和 0.25，结果为 3。权重应在完整可用的基因 TSS 目录上确定，再进行距离筛选。公式不能恢复未注释的启动子。

**证据状态。** $Q=\min\{Q_A,Q_C,Q_T,Q_{\mathrm{cat}}\}$ 分别考虑活性、接触、TSS 和候选目录质量。默认阈值为 0.5；四项均达标且基因没有未评分候选，才能标记 `sufficient_input_evidence`。分数缺失或活性/TSS 质量为零时标记 `insufficient`；其他情况为 `provisional`。未知质量保留为未知，$Q$ 不是调控概率。

**未评分与残余支持。** 若原始支持为 2、1、`NA`，且 $U$ 未知，分数为 2/3、1/3、`NA`，并报告一个未评分候选。若独立证据支持 $U=1$，则变为 1/2、1/4、`NA`。软件不自动估计 $U$；未知时仅在计算中取零，并保留 `unassigned_mass_unknown` 标记。

靶基因分配和多 TSS 思想来自 [Hecker 等（2023），generalized ABC / STARE](https://doi.org/10.1093/bioinformatics/btad062)，PACE 采用上述加权平均实现。完整定义、例子与比较边界见[英文对照指南](ABC_COMPARISON.md)。

## 参数、默认值与设置原因

| 参数或设置 | 默认值 / 输入 | 为什么需要；在哪里设置 |
| --- | --- | --- |
| 观测掩码 $m(E,i)$ | 有限测量为 1，缺失为 0 | 区分未测和实测零；由输入值自动确定 |
| 信号尺度 $a_i$ | 示例为 1；必须为正 | 分离测定单位与测定重要性；活性 JSON 的 `scale`，真实数据需预先确定 |
| 测定权重 $w_i$ | 各 1 | 没有独立依据时采用等权；JSON 的 `weight`，不代表家养动物训练所得权重 |
| 测定质量 $q(E,i)$ | 未知，已知时为 [0,1] | 让独立质量信息影响有效活性；JSON 的 `quality_column` 指定列 |
| 接触观测 $H$、期望 $D$ | 未提供时缺失 | 用观测/期望比调整先验；`contact_observed`、`contact_expected` 列，单位和分辨率需匹配 |
| 接触可靠性 $\lambda$ | 缺少必要条件时混合权重为 0 | 低可靠性时靠近距离先验；`contact_reliability` 列，不能由接触强度直接替代 |
| 接触来源 | 默认 `distance_prior` | 区分 `matched`、`surrogate`、`distance_prior`、`unknown`；`contact_source` 列 |
| TSS 权重 $\pi$ | 不同 TSS 均分 | 避免转录本计数放大；`tss_weight` 列，独立使用证据可替代均分 |
| 分配指数 $\eta$ | 1，范围 [0,1] | 控制增强子侧分配；表格命令 `--competition-power`，0 用于消融 |
| 残余支持 $U$ | 未知，计算时取零 | 记录独立估计的遗漏支持；`unassigned_mass` 列，同一基因必须一致，单位同原始支持 |
| $Q_T$、$Q_{\mathrm{cat}}$ | 未知 | 注释质量与目录覆盖需独立评估；`tss_quality`、`catalogue_quality` 列 |
| 证据阈值 | 0.5 | 输入质量的操作性规则；`--evidence-threshold`，不改变结构分数 |
| 顺式窗口 | 严格小于 5 Mb | 候选搜索起点；预测器 `--max_distance`，表格入口不自动筛选距离 |
| 距离指数 / 尺度 / 偏移 | 1.024238616787792 / 5.9594510043736655 / 5,000 bp | 定义距离先验；直接预测器提供 `--hic_gamma`、`--hic_scale`，自定义先验可走表格入口 |
| 抑制强度 $\kappa$ | 0 | 主模型不要求甲基化或抑制性测定；可选扩展仅通过核心 API |
| 输出分数阈值 | 0.02 | 示例筛选值；`pace_filter.py --threshold`，真实分析需独立校准 |
| 峰顶扩展 / 峰数上限 | 两侧各 250 bp / 150,000 | 候选构建起点；改变它们也会改变评分背景 |
| MACS2 有效基因组大小 | 通用 YAML 中为 `2.5e9` 占位值 | 按物种、组装和比对方案确定，鸡不能直接套用哺乳动物值 |

新增的是对数据缺失、质量、多启动子、靶基因分配和证据状态的显式处理。5 Mb、约 500 bp 候选宽度、峰数上限和筛选阈值属于分析起点，不是家养动物特异最优参数。RNA 可作为表达背景或后续筛选条件，不乘入主评分公式。

[完整参数表](PARAMETERS.md)区分表格 CLI、单样本 YAML 和 Snakemake 的实际控制范围。有些兼容 YAML 字段未传递到当前预测器，需按文档指定入口调整。

## Conda 安装与预先依赖

### 安装前准备

验证平台为 **Linux x86-64**。先具备 Bash、Git 和 Conda；Windows 可在 WSL2 的 Linux 环境中运行。服务器已有 Conda 时直接使用。没有 Conda 时可安装 [Miniforge](https://github.com/conda-forge/miniforge#install)，具体下载命令见[安装指南](INSTALLATION.md#before-installation)。Python 和分析包由环境文件统一安装，无需逐个预装，也不需要 GPU。

```bash
git --version
conda --version
```

### 安装主环境

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate pace
```

如果使用下载的源码 ZIP，解压后进入包含 `environment.yml` 的目录，再执行创建环境和激活命令。

| 自动安装的包 | 环境要求 | 用途 |
| --- | --- | --- |
| Python | 3.11 | 命令行程序 |
| NumPy / pandas | ≥1.26,<3 / ≥2.2,<3 | 数值计算与表格处理 |
| PyYAML | ≥6,<7 | YAML 配置 |
| SciPy | ≥1.11,<2 | 数值与分析工具 |
| matplotlib-base | ≥3.8,<4 | 无图形桌面环境下生成 QC 图 |
| pyBigWig | ≥0.3.22,<0.4 | bigWig 区间信号读取 |
| bedtools | ≥2.31,<3 | 区间操作与 reads 计数 |
| samtools | ≥1.19,<2 | BAM 与参考序列索引 |
| pytest | ≥8,<10 | 软件测试 |
| pip | ≥24 | 可选 Python 包 |

环境文件已设置 conda-forge、bioconda 和严格频道优先级。已有定量表、候选峰或处理后的信号文件时，主环境即可运行。

### 需要时安装附加功能

从已比对 reads 自动调用 MACS2、运行 Snakemake 时，使用包含全部主环境依赖的另一套环境：

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f workflow/envs/pace-env.yml
conda activate pace-workflow
snakemake --version
macs2 --version
```

它额外包含 MACS2 2.x、Snakemake 7.32.4、PuLP 2.7、setuptools <81，版本限制用于保持兼容性。两套环境按任务选择，无需同时安装。软件不执行 FASTQ 比对或原始 Hi-C 矩阵构建。

只有使用相应二进制接触文件时，才在所用环境中安装对应读取器：

```bash
# 读取 .hic：
conda install --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority hic-straw

# 读取 .cool：
conda install --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority cooler
```

BEDPE 不需要额外读取器。读取器不会自动提供接触期望或可靠性评估。

### 检查安装

```bash
python -c "import numpy, pandas, yaml, scipy, matplotlib, pyBigWig; print('Imports OK')"
bedtools --version
samtools --version
python scripts/pace.py --help
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

预期导入检查打印 `Imports OK`，48 项测试通过，小数据流程以 `PASS` 结束。确切 Linux 构建也可用随附的 [Conda 锁定文件](INSTALLATION.md#recreate-the-tested-linux-builds)安装。

## 具体使用方法

以下命令在项目根目录、`pace` 或 `pace-workflow` 环境中执行。

### 方法一：已有定量表

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/quickstart/predictions.tsv
```

输入表每行是一个增强子–基因–TSS 组合。必要列为 `chr`、`start`、`end`、`TargetGeneEnsemblID`、`TargetGeneTSS`、`contact_prior`，再加活性列或 JSON 指定的测定列。`TargetGeneEnsemblID` 是稳定基因 ID 字段名，不要求 ID 必须来自 Ensembl。

本例的 9 行输入汇总为 **6 条增强子–基因联系**；其中 2 条因活性缺失而保留为空值。控制台打印：

```text
Wrote 6 gene-level edges; all scores are uncalibrated.
```

查看分数、接触来源和证据状态：

```bash
python - <<'PYCODE'
import pandas as pd
p = pd.read_csv('results/quickstart/predictions.tsv', sep='\t')
print(p[['TargetGeneEnsemblID', 'start', 'PACE.Score',
         'contact_state', 'evidence_status']].to_string(index=False))
PYCODE
```

[可运行情景示例](WORKED_EXAMPLES.md)进一步展示缺少 H3K27ac、缺少接触和接触实测为零的不同结果。

### 方法二：已有候选峰与基因组信号

一条命令运行随附的合成 reads 示例：

```bash
bash example/run_example_direct.sh
```

输出在 `example/results/Example_Sample/`，完整表含 **12,000 条候选联系**。若要逐步替换输入，按[教程](TUTORIAL.md#step-by-step-genomic-file-example)依次执行：

| 步骤 | 程序 | 主要输入和输出 |
| --- | --- | --- |
| 1. 准备 TSS | `scripts/prepare_tss.py` | 同组装 GTF → 去重后的 TSS 表 |
| 2. 构建候选 | `workflow/scripts/pace_candidate_regions.py` | narrowPeak、染色体长度 → 候选 BED |
| 3. 定量活性 | `workflow/scripts/pace_neighborhoods.py` | 候选、注释、ATAC/DNase、可选 H3K27ac → `EnhancerList.txt`、`GeneList.txt` |
| 4. 评分 | `workflow/scripts/pace_predict.py` | 增强子与基因表、可选合格接触 → 未过滤预测表 |
| 5. 筛选与 QC | `pace_filter.py`、`pace_metrics.py` | 预测表 → 筛选表、QC 汇总和图 |

没有 H3K27ac 时可省略对应文件选项。需要明确记录“计划测定但缺失”、固定信号尺度或局部质量时，使用定量表入口。只有距离先验时仍可评分，结果通常为 `provisional`。

### 方法三：从已比对 reads 自动运行

激活 `pace-workflow`，先检查任务，再运行随附示例：

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile example/config.yaml --cores 2 --dry-run all

snakemake --snakefile workflow/Snakefile \
  --configfile example/config.yaml --cores 2 all
```

准备真实样本时，复制 `config/config.yaml` 和 `config/config_biosamples.tsv`，修改样本表路径、参考文件和输出目录。保留完整样本表表头；每个样本提供 DHS/ATAC 中的一种。详细字段及工作流限制见[Snakemake 教程](TUTORIAL.md#optional-snakemake-workflow)。

### 替换为实际物种和组织

1. 使用同一组装的 FASTA、GTF、染色体长度、峰和信号，统一 `1` / `chr1` 等名称及稳定基因 ID。
2. 在评分前固定候选集合、TSS 目录和纳入规则。重复 TSS 应去重；未注释的基因和启动子不会由评分补全。
3. 记录测定归一化、重复样本处理和信号尺度。reads/kb 不等于文库深度归一化。
4. 有 Hi-C 时提供兼容的期望接触、可靠性和来源；没有时保留距离先验。借用组织的接触要标记 `surrogate`。
5. 猪、牛、鸡分别使用对应的参考文件和有效基因组大小；鸡还需一致处理微染色体及 scaffold 的纳入范围。
6. 保留完整预测，再根据独立验证确定筛选阈值。`NA` 表示缺失，`0` 表示实测零，不能相互替代。

## 读懂结果

| 字段 | 含义 |
| --- | --- |
| `PACE.Score` | 给定候选背景下的相对支持；空值代表无法评分 |
| `contact_gene` / `target_share` / `raw_support` | 调整后的基因接触、靶分配比例、原始支持 |
| `TargetGeneTSSs` / `n_tss` | 纳入汇总的不同 TSS |
| `contact_state` | 距离先验、匹配/借用观测，或未满足质量条件的观测 |
| `evidence_status` / `evidence_reasons` | 输入证据是否充分及具体原因 |
| `score_scope` / `unscored_candidates` | 残余支持状态、当前基因内未评分候选数 |

分数不是调控概率，`sufficient_input_evidence` 也不是功能验证标签。示例阈值 0.02 不是家养动物通用的 FDR 阈值。保留未过滤预测表，再做筛选；否则可能改变归一化背景或隐藏未评分候选。

## 保存环境与引用

```bash
mkdir -p results/provenance
git rev-parse HEAD > results/provenance/pace_commit.txt
conda env export --no-builds > results/provenance/environment.yml
conda list --explicit > results/provenance/conda-explicit.txt
```

同时保存输入来源、组装和注释版本、配置、信号尺度、接触 QC 及完整结果。合成示例用于验证程序行为，生物学证据范围见 [Validation](../VALIDATION.md)。

作者与维护者：**shenlinyong，申林用（Linyong Shen），西北农林科技大学**。引用 PACE 仓库和使用的确切提交；当前仓库未指定 PACE 论文 DOI。涉及原始 ABC 和 generalized ABC 组件时同时引用上述方法论文。问题可提交至 [GitHub Issues](https://github.com/shenlinyong/PACE/issues)。
