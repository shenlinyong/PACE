# PACE

**用于家养动物增强子—基因调控研究的活性与接触整合软件。**

PACE 根据元件活性、多个启动子的综合接触和可选分配项，计算候选元件在同一基因调控背景中的相对支持。作者与维护者：**申林用（Linyong Shen，shenlinyong），西北农林科技大学**。

**第一次使用请读[完整中文说明书](docs/USER_GUIDE.zh-CN.md)**：包括安装、可直接运行的示例、真实数据准备、三种模式完整 YAML、逐项参数、结果解释和多组学接口。[English](README.md)。

## 下载安装

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
PACE --version
```

没有 Conda 时，可使用 Python 3.11 或更新版本创建环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[io,ml]'
```

Docker 用户从源码构建：

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  measured --config examples/measured/config.yaml --out results/docker_measured
```

[安装手册](docs/INSTALLATION.md)还包括 wget 下载、CPU PyTorch、Docker 序列扩展和故障排查。当前说明使用源码安装，不假设已发布 PyPI 包、公共 Docker 镜像或真实物种权重。

## 先跑通示例

```bash
PACE measured --config examples/measured/config.yaml --out results/measured
PACE hybrid --config examples/hybrid/config.yaml --out results/hybrid
PACE genome --config examples/genome_only/config.yaml --out results/genome
```

三个命令使用仓库中的合成数据，可离线运行，不需 GPU。输出目录不能已存在。真实研究需要替换数据和适用模型；示例运行成功不代表已通过家养动物生物学验证。

## 三种数据条件

| 模式 | 活性来源 | 接触来源 | 结果定位 |
|---|---|---|---|
| measured | 合格实测的固定 ATAC/DNase/H3K27ac 组合 | 实测或明确声明的适用先验 | 实测支持的调控候选 |
| hybrid | 合格实测与适用定量预测；通过校准后融合 | 实测、先验或可靠性收缩 | 混合证据预测 |
| genome_only / genome | 参考或个体序列的定量预测 | 适用距离先验 | 指定组织条件下的遗传调控潜能 |

“新个体只有 WGS”可以利用已经验证的适用模型；“整个物种没有可用功能数据”不能凭本软件自动解决。没有合适序列模型时，有实测数据就从 measured 开始。

## 按你手里的数据开始

PACE 不要求固定数量的组学层或生物学重复。先按数据条件选模式，再把每个样本的测量登记到统一表格中。

| 你已有的数据 | 先选 | 主输入 | 可选输入 |
|---|---|---|---|
| ATAC-seq、DNase-seq 或 H3K27ac，加启动子接触数据 | `measured` | `observed_activity.tsv`、`observed_contacts.tsv`、`samples.tsv` | RNA-seq、其他组蛋白、CTCF、甲基化 |
| 上述实测数据，加序列模型或个体基因组 | `hybrid` | 实测表，加 FASTA/VCF 和序列模型 | RNA-seq、其他组蛋白、CTCF、甲基化 |
| 参考/个体基因组，加已验证的序列模型 | `genome` | FASTA、VCF/BCF、callable 位点、ploidy、序列模型 | RNA-seq、其他组蛋白、CTCF、甲基化 |

主分数一次运行只使用固定的 ATAC、DNase、H3K27ac、ATAC+H3K27ac 或 DNase+H3K27ac 组合。RNA-seq、H3K4me1/3、H3K27me3、H3K9me3、CTCF 和 WGBS/RRBS 会保留为带名称和来源的注释，或进入单独验证的机器学习特征，不会被软件偷偷乘进主分数。详见[多组学接口](docs/MULTIOMICS.md)。

生物学重复可以是 1、2、3 个或更多。`samples.tsv` 为每个生物学重复和技术重复保留独立 `sample_id`，活性、接触、表达和甲基化表使用相同的 ID。软件会先在生物学重复内处理技术重复，再保留供体和 assay 信息进行评分。字段见[输入字典](docs/data_dictionary.md#samples.tsv)。

仓库里已经有可以直接运行的真实示例文件。它们使用合成样本 ID `S_ATAC`、`S_H3K27ac` 和 `S_Hi-C`，因此下面的命令可以在克隆后直接执行：

```bash
PACE measured --catalog-dir examples/measured \
  --samples examples/measured/samples.tsv \
  --activity examples/measured/observed_activity.tsv \
  --contacts examples/measured/observed_contacts.tsv \
  --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --contact-scale toy_contact \
  --no-include-promoters --out results/measured_direct
```

自己的研究数据再替换这些路径。建议文件名写出实际个体、组织、检测层和重复号，例如 `animal_001_rep2_liver_ATAC.tsv`、`animal_001_rep2_liver_HiC_5kb.tsv`；PACE 真正校验的是表头、`sample_id` 和来源元数据，不是文件名本身。

## 总公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

分子是当前元件的“活性 × 综合接触 × 可选分配”；分母是同一基因可评分候选的同类支持总和。A_star 是最终采用的活性，Cbar 是多物理 TSS 的加权接触，B 是该元件偏向目标基因的接触份额。默认实际分配指数为零，只有合格功能数据通过独立分组验证才自动采用学习结果。

[完整公式说明](docs/FORMULA.zh-CN.md)给出**全部展开的通用总公式、三种模式分别展开的总公式、每个符号和数值算例**。缺失是 NA，真实零是 0；分数是相对支持，不能解释为因果概率或基因表达贡献比例。

## 其他表观组学与结果

软件支持 RNA-seq、ATAC/DNase、H3K27ac、H3K4me1、H3K4me3、H3K27me3、H3K9me3、CTCF、WGBS/RRBS 和 Hi-C。[多组学手册](docs/MULTIOMICS.md)列出接口、角色、甲基化覆盖和启动子窗口、CTCF motif 方向及 ML 接入。附加层默认保留为可追溯注释，不随意乘入主分数。

首先查看 `gene_summary.tsv` 的候选覆盖、`qc_report.json` 的缺失原因，再看 `scores.tsv.gz` 的分数。`partial` 代表只对可评分子集归一化。`resolved_activity.tsv`、`resolved_contacts.tsv` 和证据目录说明每个值来自哪里。`run_manifest.json` 与 `resolved_config.yaml` 用于重现结果。

多物种参考版本、表观数据不齐、Hi-C 分辨率差异、多 TSS、供体重复及鸡的性染色体倍性等都有明确处理规则。这些设计使框架便于适配家养动物；具体物种、组织和品种的预测效果仍需独立实验验证。本仓库**没有附带已验证的真实家养动物序列权重**。

[完整中文说明书](docs/USER_GUIDE.zh-CN.md) · [全部参数](docs/parameters.md) ·
[输入字段](docs/data_dictionary.md) · [训练](docs/training.md) ·
[差异比较](docs/comparison.md) · [局限](docs/limitations.md)。

## 新手完整流程：从环境到结果

### 第一步：建立独立环境

推荐 Conda。不要直接使用系统 Python 或另一个项目的 NumPy。

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
python --version
python -c "import numpy, yaml; print(numpy.__version__, yaml.__version__)"
python -m pip install -e '.[io,ml]'
PACE --version
```

Python 应为 3.11 或更新版本。若 `import numpy` 报二进制错误，先运行 `conda activate pace`，再检查 `which python` 和 `which PACE` 是否都指向 `pace` 环境。没有 Conda 时，可以使用 `python3.12 -m venv .venv`、`source .venv/bin/activate` 和 `python -m pip install -e '.[io,ml]'`。

环境成功后先做一次离线自检：

```bash
PACE --help
PACE measured --help
PACE demo --regime measured --out results/check_measured
```

### 第二步：按手里的数据选模式

| 你手里的数据 | 选择 | 关键输入 | 不要误解为 |
|---|---|---|---|
| ATAC/DNase/H3K27ac 和启动子接触实测 | `measured` | 活性表、接触表、候选目录、样本表 | 不是功能验证概率 |
| 上述实测加适用序列模型或个体基因组 | `hybrid` | measured 数据 + FASTA/VCF + 序列模型 | 不是任意实测/预测平均 |
| 主要是基因组和已验证模型 | `genome` | FASTA、VCF/BCF、callable、ploidy、序列模型、接触先验 | 不是已经测得的 ATAC/RNA |

RNA-seq、H3K4me1/3、H3K27me3、H3K9me3、CTCF 和 WGBS/RRBS 都可以加入，但默认作为带来源的注释或独立 ML 特征，不会自动乘进主分数。一个项目可以有 1、2、3 个或更多生物学重复；每个重复在 `samples.tsv` 中有独立 `sample_id`，并在相应数据表中复用该 ID。

### 第三步：情况 A——只有实测组学

仓库中真实存在的 measured 文件包括 `examples/measured/samples.tsv`、`observed_activity.tsv` 和 `observed_contacts.tsv`。下面命令可以直接复制：

```bash
PACE measured \
  --catalog-dir examples/measured \
  --samples examples/measured/samples.tsv \
  --activity examples/measured/observed_activity.tsv \
  --contacts examples/measured/observed_contacts.tsv \
  --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --contact-scale toy_contact \
  --no-include-promoters --out results/measured_direct
```

真实 measured 项目至少需要：

```text
samples.tsv              每个生物学/技术重复和供体
observed_activity.tsv    ATAC、DNase、H3K27ac 的定量窗口
observed_contacts.tsv    Hi-C/Prom-Hi-C 的元件—启动子接触
units/promoters/candidates.tsv  固定候选目录
sources.tsv/evidence.tsv         来源和证据链
```

### 第四步：情况 B——实测加序列预测

hybrid 在 measured 文件基础上增加参考基因组、个体变异、可调用区域、倍性和序列模型；只有有匹配校准器时才做实测—预测融合。

```bash
PACE hybrid \
  --catalog-dir examples/hybrid \
  --samples examples/hybrid/samples.tsv \
  --activity examples/hybrid/observed_activity.tsv \
  --contacts examples/hybrid/observed_contacts.tsv \
  --reference examples/hybrid/genome.fa --vcf examples/hybrid/sample.vcf \
  --callable examples/hybrid/callable.bed --ploidy examples/hybrid/ploidy.tsv \
  --sequence-model examples/hybrid/models/sequence \
  --fusion-model examples/hybrid/models/fusion \
  --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --contact-scale toy_contact \
  --no-include-promoters --sample-id toy_animal --individual-id toy_animal \
  --out results/hybrid_direct
```

实测样本和 VCF 个体必须对应同一个研究对象，物种、assembly、组织、窗口和单位必须与模型清单一致。

### 第五步：情况 C——主要只有基因组

genome 不要求 ATAC 或 RNA 文件，但不能只给 FASTA；还需要适用序列模型和接触先验。个体分析还需要 VCF/BCF、callable 和 ploidy。

```bash
PACE genome \
  --catalog-dir examples/genome_only \
  --reference examples/genome_only/genome.fa --vcf examples/genome_only/sample.vcf \
  --callable examples/genome_only/callable.bed --ploidy examples/genome_only/ploidy.tsv \
  --sequence-model examples/genome_only/models/sequence \
  --contact-prior examples/genome_only/models/contact \
  --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --contact-scale toy_contact \
  --no-include-promoters --sample-id toy_animal --individual-id toy_animal \
  --out results/genome_direct
```

### 第六步：检查输出

每次运行先看 `qc_report.json` 和 `gene_summary.tsv`，确认实际有多少候选进入分母；再看 `scores.tsv.gz` 的 `pace_score`、`A_used`、`Cbar`、`support` 和 `reason`。`resolved_activity.tsv`、`resolved_contacts.tsv` 说明每个值来自实测、预测、融合还是先验；`run_manifest.json` 记录配置、输入和模型身份。

## 真实数据到底要整理成什么样

PACE 从已经处理好的定量表开始，不负责 FASTQ 比对、peak calling、变异检测、相位推断或通用 liftover。先把不同来源整理成下面的规范表，再运行评分。

### `samples.tsv`：样本、供体和重复

至少包含以下字段：

```text
sample_id  donor_id  assay  biological_replicate  technical_replicate  species  assembly  context_id  source_id
```

例如同一只动物的三个肝脏生物学重复可以写成：

```text
animal_001_liver_rep1_ATAC  animal_001  ATAC     1  1  cattle  ARS-UCD1.2  liver  lab_batch_01
animal_001_liver_rep2_ATAC  animal_001  ATAC     2  1  cattle  ARS-UCD1.2  liver  lab_batch_02
animal_001_liver_rep3_ATAC  animal_001  ATAC     3  1  cattle  ARS-UCD1.2  liver  lab_batch_03
animal_001_liver_rep1_HiC   animal_001  Hi-C     1  1  cattle  ARS-UCD1.2  liver  lab_batch_01
```

不要求每种 assay 都有相同的重复数，但缺少的样本必须在质量状态和 QC 中保持可见。技术重复应使用同一个生物学重复号和不同的 `technical_replicate`，不要把技术重复伪装成更多动物。

### `observed_activity.tsv`：主活性层

每一行是一个固定元件、一个样本和一个 assay：

```text
element_id  sample_id  assay  signal  measurement_status  callable_fraction  unit  normalization_id  window_id
E0001        animal_001_liver_rep1_ATAC  ATAC     12.4  observed  0.98  CPM  atac_norm_v1  grid:500:mean
E0001        animal_001_liver_rep1_H3K27ac H3K27ac 8.1 observed 0.97  CPM  chip_norm_v1  grid:500:mean
```

主活性组合只能从 ATAC、DNase、H3K27ac 中选择单层或支持的双层组合，并在一次运行中保持不变。RNA-seq 不放在这里；它进入 `expression.tsv`。H3K4me1、H3K4me3、H3K27me3、H3K9me3 和 CTCF 通常进入 `features.tsv` 或相应准备接口。

### `observed_contacts.tsv`：元件—启动子接触

每一行对应一个元件、一个物理启动子和一个接触样本：

```text
element_id  promoter_id  sample_id  contact_value  measurement_status  bin_pair_id  scale  resolution  source_id
E0001        ENSG000001_P1  animal_001_liver_rep1_HiC  0.34  observed  bin_102:bin_205  balanced_HiC  5000  hic_batch_01
```

`resolution`、`scale`、`normalization_id` 和 `balancing` 必须和先验或其它接触样本兼容。p 值、相关系数和未声明背景的 log(O/E) 不能直接当作接触值。

### RNA-seq、甲基化和其它组学

| 数据 | 文件 | 最小身份字段 | 在主流程中的默认角色 |
|---|---|---|---|
| RNA-seq | `expression.tsv` | `gene_id`, `sample_id`, `tpm`, `status` | 基因表达注释或独立 ML 特征 |
| WGBS/RRBS | `methylation.tsv` | `chrom`, `dyad_start0`, `methylated_count`, `total_count`, `sample_id`, `assay` | CpG 覆盖和甲基化状态 |
| 组蛋白/CTCF | `features.tsv` 或 bigWig/BED 准备输出 | `entity_id`, `feature_name`, `value`, `evidence_id`, `status` | 元件、启动子或结构注释 |

这些数据可以有不同的重复数，但每一行都要保留 `sample_id` 和来源。PACE 不会因为用户提供了 RNA-seq 就自动把 TPM 乘到主分数上；这避免了把表达注释误读为因果效应。

## 从原始处理结果到规范表

真实轨道先用准备命令转换，不能把 bigWig、BAM 或 xlsx 直接改名为 TSV：

```yaml
kind: bigwig
track: data/animal_001_rep1_liver_ATAC.bw
units: prepared/catalog/units.tsv
sample_id: animal_001_liver_rep1_ATAC
assay: ATAC
unit: normalized_signal
normalization_id: atac_norm_v1
window_id: grid:500:mean
missing_is_measured_zero: false
minimum_callable_fraction: 0.8
```

```bash
PACE prepare --config prepare_atac.yaml --out prepared/animal_001_rep1_atac
```

对第二、第三个重复重复准备步骤，最后把输出行合并到同一张 `observed_activity.tsv`，保留一次表头和每行的 `sample_id`。Hi-C/cooler、RNA、甲基化和 CTCF 的准备方式见[输入准备](docs/input_preparation.md)和[多组学接口](docs/MULTIOMICS.md)。

## 常见问题先看这里

| 现象 | 先检查 |
|---|---|
| `numpy.core._multiarray_umath` | 是否激活了 `pace` 环境；`which python` 和 `which PACE` 是否来自同一环境 |
| `Output already exists` | 结果目录是保护性的，换一个新目录，不要删除旧结果 |
| `Contact scale differs` | 实测接触、先验、配置中的 `contact.scale` 是否相同 |
| `missing required column` | 查看[数据字典](docs/data_dictionary.md)，不要用文件名代替字段 |
| 活性变成 NA | 检查 `measurement_status`、`callable_fraction`、单位和 `activity.panel` 是否匹配 |
| 只有 WGS 但没有模型 | genome 模式不能从 FASTA 自动训练组织特异模型，需要适用权重和接触先验 |

公式中的每个因子、零值与缺失规则以及三种模式的展开式见[中文公式说明](docs/FORMULA.zh-CN.md)。
