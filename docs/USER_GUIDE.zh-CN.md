# PACE 完整使用说明

PACE 用于整合增强子活性和增强子—启动子接触，计算候选调控元件对同一个基因的相对支持。本文按“安装 → 跑通示例 → 准备自己的数据 → 选择模式 → 检查结果”的顺序介绍。作者与维护者：申林用（Linyong Shen，shenlinyong），西北农林科技大学。

这里的“支持”是候选元件之间的相对份额。PACE 分数高，不等于该联系已经得到功能实验验证，也不等于该元件贡献了相同比例的基因表达。软件会保留原始支持、证据来源、缺失原因和实际分母，供研究者判断。

## 1. 第一次使用：下载、安装和运行

### 1.1 从 GitHub 下载

Linux 服务器执行：

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
git rev-parse HEAD
```

保存最后一条命令输出的提交号，论文分析应固定这一版本。没有 Git 时可下载源码压缩包：

```bash
wget -O PACE-main.tar.gz https://github.com/shenlinyong/PACE/archive/refs/heads/main.tar.gz
tar -xzf PACE-main.tar.gz
cd PACE-main
```

### 1.2 推荐的 Conda 安装

```bash
conda env create -f environment.yml
conda activate pace
PACE --version
PACE measured --help
```

需事先安装 Conda。该环境包含核心程序及常用基因组文件读取依赖。真实序列 CNN 需要另外安装 PyTorch 和序列扩展；[安装手册](INSTALLATION.md)给出了 CPU 命令、Python venv、独立安装目录和 Docker 三种替代方式。

Docker 从源码本地构建，容器入口已经是 PACE：

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  measured --config examples/measured/config.yaml --out results/docker_measured
```

### 1.3 先运行三个自带示例

在仓库根目录执行：

```bash
PACE measured --config examples/measured/config.yaml --out results/demo_measured
PACE hybrid --config examples/hybrid/config.yaml --out results/demo_hybrid
PACE genome --config examples/genome_only/config.yaml --out results/demo_genome
```

这些命令可以直接运行，数据和示例权重已在仓库中，不需要另找文件，也不需要 GPU。它们全部是**合成数据**，用于熟悉命令和检查安装。输出目录必须尚不存在；重复运行时换一个名称。

真实研究请使用自己的数据和适用模型，并把 `execution_profile` 设为 `research`。把合成模型的 `is_synthetic` 改为 false 不能让它成为真实模型。仓库没有附带已验证的鸡、猪、牛等物种权重。

## 2. 应该选择哪一种模式

| 情况 | 模式 | 活性从哪里来 | 接触从哪里来 | 可以解释为什么 |
|---|---|---|---|---|
| 有质量合格的 ATAC/DNase/H3K27ac 等测量 | `measured` | 所选固定活性组合的实测信号 | 实测 Hi-C；也可显式配置适用先验 | 实测证据支持的调控候选排序 |
| 部分测量缺失或不足，并有适用序列模型 | `hybrid` | 合格实测与定量序列预测；有匹配校准器才融合 | 实测、先验或二者按已声明可靠性组合 | 混合证据下的候选排序 |
| 新个体只有基因组数据，但已有适用功能模型 | `genome_only`，命令别名 `genome` | 该个体序列的定量预测 | 已拟合的距离接触先验 | 指定组织背景下的遗传调控潜能 |

三种模式共享同一个评分公式。模式差别是证据来源，不能把“实测”和“预测”混成同一种实验结果。`measured` 允许 H3K27ac 单层运行；没有 ATAC 时可以明确选择这个组合，不需要伪造 ATAC 数值。

**没有合适序列模型时：**有实测数据就使用 measured；只有 WGS 则不能从这个仓库直接得到可靠的组织特异调控预测。某个新个体没有表观数据，与整个物种都没有训练和验证数据，是两种不同情况。

### 2.1 一份、两份还是三份重复都可以

PACE 不把重复数写死。每个生物学重复、技术重复和供体都在 `samples.tsv` 中占一行，并在其它表中通过相同的 `sample_id` 关联：

```text
sample_id        donor_id   assay       biological_replicate   technical_replicate
animal_001_rep1  animal_001 ATAC        1                      1
animal_001_rep2  animal_001 ATAC        2                      1
animal_001_rep3  animal_001 H3K27ac     3                      1
```

`observed_activity.tsv` 可以同时包含这些样本的 ATAC、DNase 或 H3K27ac 行，`observed_contacts.tsv` 可以同时包含对应的 Hi-C/Prom-Hi-C 样本行。RNA-seq、H3K4me1/3、H3K27me3、H3K9me3、CTCF 和 WGBS/RRBS 也保留自己的 `sample_id` 与 `assay`。软件会先在生物学重复内处理技术重复，再按照目标层级保留供体信息；不会因为恰好有三个重复就自动声称统计功效足够。

### 2.2 其它组学放在哪里

ATAC/DNase/H3K27ac 是主活性组合的候选层；Hi-C/Prom-Hi-C 是主接触证据。RNA-seq 进入 `expression.tsv`，甲基化进入 `methylation.tsv`，其它组蛋白和 CTCF 进入 `features.tsv` 或相应的准备接口。它们默认用于注释和可追溯性，只有在独立功能标签、明确特征契约和验证通过时才进入 ML 分析，不会被自动乘进 PACE 主公式。字段和示例见[多组学说明](MULTIOMICS.md)。

## 3. 总公式与阅读顺序

```math
\mathrm{PACE}(E,G)=\frac{
A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}.
```

按三个步骤理解：先确定每个元件的活性，再确定它对该基因的综合接触，最后把它的支持除以该基因可评分候选的支持总和。

| 符号 | 通俗解释 |
|---|---|
| E、G | 当前调控元件和目标基因 |
| e | 分母中逐个遍历的候选元件 |
| A_star | 最终用于计算的活性，可能来自实测、预测或已校准融合 |
| Cbar | 把同一个基因不同物理 TSS 的接触按固定权重汇总 |
| B | 该元件对该基因的接触占其所有候选目标基因接触的份额 |
| eta_used | 实际使用的可选分配参数；默认数值为 0，有合格功能证据才自动学习并验证 |
| E_score | 在预先确定的候选背景中，本次确实可计算支持的元件 |

以 ATAC+H3K27ac 为例，活性是两种合格信号的几何平均：

```math
A_\star(E)=\sqrt{x_{\star,\mathrm{ATAC}}(E)\,x_{\star,\mathrm{H3K27ac}}(E)}.
```

两个转录本共用同一个 TSS，只计一个物理启动子。有多个 TSS 时，用预先提供的启动子权重；没有可信使用比例时可以明确设为等权。

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)
\{r(E,t)C_{\mathrm{obs}}(E,t)+[1-r(E,t)]C_{\mathrm{prior}}(E,t)\}.
```

这里 pi 是同一基因的 TSS 权重，r 是接触中实测来源的权重；二者含义不同。r=1 只用实测，r=0 只用先验。程序不会仅凭测序文件的存在自动认定数据可信。

[完整数学说明](FORMULA.zh-CN.md)依次列出了：简写总公式、全部展开的通用总公式、measured 详细总公式、hybrid 详细总公式、genome_only 详细总公式、逐个符号及零值/缺失规则。hybrid 的实测—预测权重在单个检测层的 log1p 空间校准；genome_only 先平均各条染色体拷贝的预测信号，再计算活性。两者都不是把几个最终 PACE 分数相加。

基因表达权重没有再乘入主分数。同一基因的 TPM 权重同时乘到分子和分母会抵消；只在归一化后乘则改变分数含义。因此 RNA 作为独立注释或经过验证的机器学习特征输出。

## 4. 真实数据运行前需要哪些文件

软件从**已经处理的基因组数据**开始，不执行 FASTQ 比对、峰识别、变异检测、相位推断或通用坐标转换。所有文件必须使用同一物种、参考组装版本和染色体命名。

| 文件 | 是否需要 | 内容 |
|---|---|---|
| `units.tsv` | 所有模式必需 | 非重叠固定评分单元，默认 500 bp |
| `promoters.tsv` | 所有模式必需 | 基因、物理 TSS、方向、TSS 权重 |
| `candidates.tsv` | 所有模式必需 | 预先确定的元件—基因候选集合 |
| `samples.tsv` | 使用真实测量时必需 | 样本、个体、检测类型、生物/技术重复与背景 |
| `sources.tsv` | 有来源引用时必需 | 文件或 accession、处理方法、归一化、校验和 |
| `evidence.tsv` | 输入引用证据或导入解析结果时必需 | 证据身份和来源链；不能随意编造测量样本 |
| `observed_activity.tsv` | measured；hybrid 可用 | 按固定单元汇总的定量活性及测量状态 |
| `observed_contacts.tsv` | 采用实测接触时 | 元件到各 TSS 的接触值、分辨率和尺度 |
| 序列模型目录 | 需要定量序列预测时 | manifest 与真实模型权重；必须匹配物种、组织、输出量纲 |
| 接触先验目录 | 使用先验或收缩时 | 同尺度、同分辨率且适用背景的已拟合先验 |
| `expression.tsv`、甲基化、`features.tsv` | 可选 | RNA 和其他表观组学证据 |

各表必须是带表头的 UTF-8 TSV；不要把 Excel xlsx 直接重命名为 tsv。`NA` 表示缺失，0 表示真实测得或预测的零。BED 和内部区间使用 0-based、左闭右开坐标；`tss0` 是 0-based 位点。[完整字段字典](data_dictionary.md)给出所有必需表头。

建议项目中分开存放：`data/` 原始处理后文件，`prepared/` 规范表，`models/` 模型资产，`configs/` 配置，`results/` 输出。以下真实数据配置假设 YAML 放在项目根目录，因此路径直接写 `prepared/...`。如果移动到 configs 目录，路径需相应加 `../`。

## 5. 从 BED/GTF、bigWig 和 Hi-C 准备规范输入

### 5.1 建立候选目录

把下面保存为 `prepare_catalog.yaml`，替换为自己的文件名：

```yaml
kind: catalog
bed: data/enhancer_atlas.bed
gtf: data/annotation.gtf
chrom_sizes: data/chrom_sizes.tsv
source_id: chicken_liver_atlas
width: 500
offset: 0
radius: 5000000
include_promoters: true
```

```bash
PACE prepare --config prepare_catalog.yaml --out prepared/catalog
```

`chrom_sizes.tsv` 有 `chrom`、`length` 两列。GTF 需要 transcript 记录。`radius` 是候选搜索距离，不是鸡的已验证最优距离，应在研究设计中确定并检查敏感性。启动子单元默认进入候选背景。参考序列边缘放不下完整窗口的单元会被报告。

### 5.2 提取活性轨道

`prepare_atac.yaml`：

```yaml
kind: bigwig
track: data/animal1_ATAC.bw
units: prepared/catalog/units.tsv
sample_id: animal1_ATAC
assay: ATAC
unit: normalized_signal
normalization_id: frozen_signal_protocol
window_id: grid:500:mean
missing_is_measured_zero: false
minimum_callable_fraction: 0.8
```

```bash
PACE prepare --config prepare_atac.yaml --out prepared/animal1_atac
```

H3K27ac 使用同一命令，修改轨道、sample_id、assay 和输出目录。0.8 只是这里演示的阈值，应按实际轨道协议设置。bigWig 未存储的区间不一定是真实零，不能不加判断就打开 `missing_is_measured_zero`。

把各样本输出的 `observed_activity.tsv` 合并成一张保留一次表头的 TSV。使用 Python 标准库即可：

```bash
python - <<'PY'
import csv
from pathlib import Path
paths = [Path('prepared/animal1_atac/observed_activity.tsv'),
         Path('prepared/animal1_h3k27ac/observed_activity.tsv')]
with Path('prepared/observed_activity.tsv').open('w', newline='') as out:
    writer = None
    for path in paths:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if writer is None:
                writer = csv.DictWriter(out, reader.fieldnames, delimiter='\t')
                writer.writeheader()
            elif reader.fieldnames != writer.fieldnames:
                raise ValueError(f'Headers differ: {path}')
            writer.writerows(reader)
PY
```

需先对 H3K27ac 执行准备命令，生成第二个路径。不要重复输入同一 sample/element/assay。各检测层可以有自己的量纲，但同一层在所有元件、样本和所用模型之间必须采用同一个已声明归一化协议。

### 5.3 建立 E–TSS 查询对并提取 Hi-C

`element_promoter_pairs.tsv` 每行包含 `element_id,promoter_id,chrom,anchor0,tss0`。从候选表与启动子表连接生成：

```bash
python - <<'PY'
import csv
from collections import defaultdict
from pathlib import Path
root = Path('prepared/catalog')
def read(name):
    with (root / name).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
units = {r['element_id']: r for r in read('units.tsv')}
promoters = defaultdict(list)
for r in read('promoters.tsv'):
    promoters[r['gene_id']].append(r)
rows = {}
for edge in read('candidates.tsv'):
    u = units[edge['element_id']]
    for p in promoters[edge['gene_id']]:
        key = (u['element_id'], p['promoter_id'])
        rows[key] = dict(element_id=u['element_id'], promoter_id=p['promoter_id'],
                         chrom=u['chrom'], anchor0=u['anchor0'], tss0=p['tss0'])
with Path('prepared/element_promoter_pairs.tsv').open('w', newline='') as out:
    writer = csv.DictWriter(out, ['element_id','promoter_id','chrom','anchor0','tss0'], delimiter='\t')
    writer.writeheader()
    writer.writerows(rows.values())
PY
```

`prepare_hic.yaml`：

```yaml
kind: cooler
contact: data/animal1.mcool::/resolutions/5000
pairs: prepared/element_promoter_pairs.tsv
resolution: 5000
balanced: true
missing_pixels_are_zero: true
sample_id: animal1_HiC
source_id: animal1_HiC_source
scale: balanced_contact_protocol_1
normalization_id: hic_norm_protocol_1
```

```bash
PACE prepare --config prepare_hic.yaml --out prepared/animal1_hic
```

`balanced` 必须符合文件实际处理方式；未存储 pixel 是否代表零也需核实。不同分辨率、raw counts、不同归一化不能直接平均。O/E 已去掉距离背景，需要先用匹配的 expected 还原为接触量；p 值不能作 contact。`.hic` 需在上游转换为可核验的 cool/mcool。

## 6. 模式一：实测数据 measured

以下是**真实项目配置模板**，需先准备所列文件；它不是仓库中开箱即用的示例。以一个鸡个体肝脏为例，保存为 `measured.yaml`：

```yaml
schema_version: pace-1
run_id: chicken_liver_animal1
regime: measured
execution_profile: research
estimand: bulk_proxy
target_level: individual
context:
  species: chicken
  assembly: GRCg7w
  context_id: liver
inputs:
  units: prepared/catalog/units.tsv
  promoters: prepared/catalog/promoters.tsv
  candidates: prepared/catalog/candidates.tsv
  region_membership: prepared/catalog/region_membership.tsv
  samples: prepared/samples.tsv
  sources: prepared/sources.tsv
  observed_activity: prepared/observed_activity.tsv
  observed_contacts: prepared/animal1_hic/observed_contacts.tsv
catalog:
  profile: canonical_grid
  width_bp: 500
  offset_bp: 0
  include_promoter_units: true
  chrom_sizes_path: prepared/catalog/chrom_sizes.tsv
  candidate_radius_bp: 5000000
activity:
  panel: [ATAC, H3K27ac]
  minimum_callable_fraction: 0.8
contact:
  mode: observed
  scale: balanced_contact_protocol_1
  resolution: 5000
  normalization_id: hic_norm_protocol_1
  balancing: balanced
  window_id: bin_pair
  near_diagonal_policy: prior_or_unresolved
  allow_prior_fallback: false
promoters:
  weights: provided
allocation:
  eta: auto
multiomics:
  mode: annotate
seed: 17
```

```bash
PACE validate --config measured.yaml
PACE measured --config measured.yaml --out results/animal1_measured
```

| 参数 | 如何选择 |
|---|---|
| `context` | 与参考序列、注释和所有来源表一致；`chicken`、`liver` 是一致性标识，不能混用不同拼法 |
| `target_level` | 单一个体用 individual；多个个体的等权平均用 population_mean，并清楚解释群体背景 |
| `activity.panel` | ATAC、DNase、H3K27ac 三种单层，或可及性+H3K27ac 双层；不支持把任意组蛋白直接塞入主活性组合 |
| `minimum_callable_fraction` | 低于阈值的观测不可用；由数据生产协议确定 |
| `contact.scale` | 与接触表一致的测量尺度标识；更名不能实现尺度转换 |
| `near_diagonal_policy` | 同 bin 等近距离接触使用匹配先验，或返回 NA；无先验时保持缺失 |
| `promoters.weights` | provided 使用可信 pi；equal 对去重后的物理 TSS 等权 |
| `allocation.eta` | 通常保持 auto，无合格校准证据时实际为零 |

如果只有 H3K27ac，把 `panel` 明确改为 `[H3K27ac]`，其余处理不变。为了使 C1/C2 生物学重复可比较，建议每个个体分别评分，检查重复一致性，再决定是否另做 population_mean 图谱；技术重复不能充当独立个体。

## 7. 模式二：混合证据 hybrid

保存为 `hybrid.yaml`。此模板同时演示序列活性与接触收缩；模型和校准器需事先训练或取得适用资产：

```yaml
schema_version: pace-1
run_id: chicken_liver_hybrid_animal1
regime: hybrid
execution_profile: research
estimand: bulk_proxy
target_level: individual
context:
  species: chicken
  assembly: GRCg7w
  context_id: liver
inputs:
  units: prepared/catalog/units.tsv
  promoters: prepared/catalog/promoters.tsv
  candidates: prepared/catalog/candidates.tsv
  samples: prepared/samples.tsv
  sources: prepared/sources.tsv
  observed_activity: prepared/observed_activity.tsv
  observed_contacts: prepared/animal1_hic/observed_contacts.tsv
catalog:
  profile: canonical_grid
  width_bp: 500
  include_promoter_units: true
  chrom_sizes_path: prepared/catalog/chrom_sizes.tsv
  candidate_radius_bp: 5000000
activity:
  panel: [ATAC, H3K27ac]
  minimum_callable_fraction: 0.8
contact:
  mode: shrinkage
  scale: balanced_contact_protocol_1
  resolution: 5000
  normalization_id: hic_norm_protocol_1
  balancing: balanced
  window_id: bin_pair
  prior_path: models/chicken_liver_contact
  reliability: 0.7
  reliability_source: held_out_contact_reliability_protocol
  allow_prior_fallback: true
sequence:
  model_path: models/chicken_liver_sequence
  max_n_fraction: 0.05
fusion:
  calibrator_path: models/chicken_liver_fusion
  quality_stratum: default
genome:
  individual_id: animal1
  sample_id: animal1
  reference_path: data/reference.fa
  variant_path: data/animal1.vcf.gz
  callability_path: data/animal1_callable.bed
  ploidy_path: data/animal1_ploidy.tsv
  unrecorded_site_policy: require_callable
allocation:
  eta: auto
multiomics:
  mode: annotate
seed: 17
```

```bash
PACE capabilities --config hybrid.yaml
PACE validate --config hybrid.yaml
PACE hybrid --config hybrid.yaml --out results/animal1_hybrid
```

| 新增参数 | 作用与要求 |
|---|---|
| `sequence.model_path` | 定量活性模型目录，物种/组装/组织/窗口/信号量纲必须匹配 |
| `fusion.calibrator_path` | 实测与预测的融合校准资产；缺少时按来源优先规则选择，不猜测融合权重 |
| `fusion.quality_stratum` | 选用校准器中实际存在的质量分层 |
| `contact.mode: shrinkage` | 用已声明可靠性组合实测接触和先验 |
| `contact.reliability` | 此处 0.7 仅演示语法，不能照抄为默认生物参数；由独立评估确定 |
| `contact.prior_path` | 接触先验，须与实测的尺度、分辨率及背景兼容 |
| `genome.*` | 个体序列重建所需参考、基因型、可调用区和倍性；实测 donor_id 需与 individual_id 一致 |

如果接触数据充分，可以改用 `contact.mode: observed`，去掉 reliability、reliability_source、prior_path 并关闭 prior fallback。如果只做参考序列预测，可省略个体 VCF、callability、ploidy、individual_id 和 sample_id，但输出不再代表该个体的基因型效应。详细训练命令见[模型训练](training.md)。

## 8. 模式三：只有新个体基因组数据 genome_only

保存为 `genome.yaml`：

```yaml
schema_version: pace-1
run_id: chicken_liver_genome_animal1
regime: genome_only
execution_profile: research
estimand: bulk_proxy
target_level: individual
context:
  species: chicken
  assembly: GRCg7w
  context_id: liver
inputs:
  units: prepared/catalog/units.tsv
  promoters: prepared/catalog/promoters.tsv
  candidates: prepared/catalog/candidates.tsv
catalog:
  profile: canonical_grid
  width_bp: 500
  include_promoter_units: true
  chrom_sizes_path: prepared/catalog/chrom_sizes.tsv
  candidate_radius_bp: 5000000
activity:
  panel: [ATAC, H3K27ac]
contact:
  mode: prior_only
  scale: balanced_contact_protocol_1
  resolution: 5000
  normalization_id: hic_norm_protocol_1
  balancing: balanced
  window_id: bin_pair
  prior_path: models/chicken_liver_contact
  allow_prior_fallback: true
sequence:
  model_path: models/chicken_liver_sequence
  max_n_fraction: 0.05
genome:
  individual_id: animal1
  sample_id: animal1
  reference_path: data/reference.fa
  variant_path: data/animal1.vcf.gz
  callability_path: data/animal1_callable.bed
  ploidy_path: data/animal1_ploidy.tsv
  phase_policy: require_phase_or_single_variant_scenario
  unrecorded_site_policy: require_callable
  sv_assessed: false
allocation:
  eta: auto
multiomics:
  mode: annotate
seed: 17
```

```bash
PACE capabilities --config genome.yaml
PACE validate --config genome.yaml
PACE genome --config genome.yaml --out results/animal1_genome
```

| 参数或输入 | 为什么需要 |
|---|---|
| 候选目录 | WGS 本身不会自动告诉软件哪些区域是组织增强子；需已有 atlas/可信候选来源 |
| 参考 FASTA | 无压缩的 A/C/G/T/N 序列，与注释和 VCF 一致 |
| VCF/BCF | 该个体规范化的变异；sample_id 必须是实际样本列名 |
| callable BED | 说明哪些没有 VCF 记录的区域可以视为参考基因型；“没有记录”不等于“没有变异” |
| `ploidy.tsv` | 有 chrom、ploidy 两列，显式说明各染色体 1 或 2 倍；鸡的性染色体不能全部默认二倍体 |
| `phase_policy` | 个体单倍型需要可信相位；无法连通的相位区块或未定相杂合变异不会被任意拼接 |
| `max_n_fraction` | 限制输入序列未知碱基比例；默认 0.05 是工程配置，需要按模型适用范围审查 |
| `sv_assessed` | 记录是否评估结构变异，不能通过设为 true 代替真实 SV 分析 |

`genome_only` 不能输入实测 activity/contact 表；若需要混合它们，应选择 hybrid。接触目前是距离先验，不是实测个体三维基因组。只有适用模型、基因型和结构检查都通过，结果才具备对应解释范围。

常规 SNV 和保持固定输出目标的短 indel 可进入序列重建。改变目标窗口、已知结构变异、BND 断点等超出支持范围时保留明确不可用状态；软件不重建复杂 CNV 剂量和全基因组重排。未被个体 GT 选中的 ALT 不应影响该个体。

需要分步重建并预测时：

```bash
PACE prepare-genome --config genome.yaml --out prepared/animal1_windows
PACE predict-sequence --config genome.yaml --out prepared/animal1_predictions
```

导入外部个体预测仍需保留原参考、VCF、callability、ploidy 和相关配置，并使用匹配的 `genome_binding_id`。该标识将预测绑定到实际准备的个体输入；软件还会重新检查结构和窗口，不能用预计算数值覆盖失败状态。它是来源一致性记录，不是外部模型准确性的证明。

## 9. RNA、甲基化及其他表观数据怎样使用

| 数据 | 主公式中的角色 | 软件入口 |
|---|---|---|
| ATAC-seq / DNase-seq | 固定活性组合中的可及性信号 | bigWig prepare → observed_activity |
| H3K27ac | 固定活性组合中的乙酰化信号 | bigWig prepare → observed_activity |
| Hi-C | 元件到物理 TSS 的接触 | cooler prepare → observed_contacts |
| RNA-seq | 基因表达注释、显式选择的 ML 特征 | expression 或 rna prepare |
| H3K4me1 | 元件状态注释，可用于独立 ML | observed_activity 中的附加层或 features |
| H3K4me3 | 启动子状态注释，可用于独立 ML | 启动子窗口量化后 features |
| H3K27me3 / H3K9me3 | 抑制性染色质注释，可用于独立 ML | bigWig 或 features；没有固定惩罚倍数 |
| CTCF ChIP-seq / CUT&Tag | 结构相关占据与 motif 方向注释 | bigWig、BED features 或外部 features |
| DNA methylation：WGBS/RRBS | 区域/启动子甲基化及覆盖状态注释 | CpG counts、methylation prepare、features |

只把文件路径写进配置，不能自动建立缺失的生物学关系。例如 CTCF 强不必然代表接触强，所有区域的甲基化也不能一律解释为抑制。完整接口、甲基化位点与 pooled 均值差别、RNA 状态和证据链见[多组学操作手册](MULTIOMICS.md)。

有独立功能标签时可训练额外分类器，输出 `pace_ml_score`，在适用校准通过时输出 `pace_ml_probability`。这些列独立于主 `pace_score`，不与主分数简单平均。模型适用范围、缺失处理、训练/校准/测试分离必须在报告中保留。

## 10. 如何看结果

```bash
python - <<'PY'
import csv, gzip, json
from itertools import islice
from pathlib import Path
out = Path('results/demo_measured')
with gzip.open(out / 'scores.tsv.gz', 'rt') as handle:
    for row in islice(csv.DictReader(handle, delimiter='\t'), 5):
        print({k: row.get(k) for k in
               ['element_id','gene_id','pace_score','A_used','Cbar','normalization_status','reason']})
with (out / 'qc_report.json').open() as handle:
    print(json.load(handle))
PY
```

| 输出文件 | 建议检查内容 |
|---|---|
| `scores.tsv.gz` | 每条边的 pace_score、A_used、Cbar、B、support、log_support、缺失原因 |
| `gene_summary.tsv` | 每个基因计划候选数、可评分数、覆盖率和分母 |
| `resolved_activity.tsv` | 每层活动最终采用实测、预测还是融合，及使用量纲 |
| `resolved_contacts.tsv` | 接触来源、prior、reliability 和无效原因 |
| `multiomics_features.tsv.gz` | 额外组学特征、状态、角色和 evidence_id |
| `evidence.tsv`、`sources.tsv` | 原始与派生证据的可追溯目录 |
| `eta_calibration.json` | 实际参数、自动回退或通过验证的依据 |
| `qc_report.json`、`report.md` | 总体覆盖、失败原因、能力边界 |
| `resolved_config.yaml`、`run_manifest.json` | 实际参数、输入/模型哈希、环境和软件身份 |

先检查覆盖与缺失，再看分数。`complete` 只表示计划候选可评分，不表示发现了全部生物学增强子；`partial` 表示分数只针对可评分子集。只有一个有效候选时它可能得到 1，但并没有因此成为已验证调控元件。所有支持为零时 PACE 是 NA，不能人为加常数凑出分数。

## 11. 比较两只动物或两个处理

两个结果必须使用相同候选背景和兼容的测量、模型与参数。`compare` 会在共同可评分元件上重新归一化，而不是直接连接两张分数表后相减。

`compare.yaml`：

```yaml
left: results/animal1_measured
right: results/animal2_measured
minimum_common_units: 2
allow_eta_difference: false
allow_evidence_difference: false
```

```bash
PACE compare --config compare.yaml --out results/animal2_minus_animal1
```

差值方向是右减左。完整差异与条件差异分开报告；分辨率、归一化或模型策略不兼容时不能当作完整个体生物效应。增强子自身支持不变，但其他元件改变，也会使它的 PACE 份额改变，因此应同时查看活性、support、基因总支持和 PACE 差值。

## 12. 为什么这个框架适配家养动物

| 家养动物分析中的实际问题 | PACE 的处理 | 仍需研究者解决的部分 |
|---|---|---|
| 物种、品种和参考版本多 | 显式绑定物种、组装、组织和模型范围 | 参考质量和跨版本坐标转换 |
| 表观检测不齐全 | 固定单层/双层活性组合；合格时使用混合证据 | 没有合适模型不能补成可信实测 |
| Hi-C 深度、分辨率和组织覆盖有限 | 接触测量合同检查、显式先验或可靠性收缩 | 适用先验拟合和独立评估 |
| 转录本/TSS 注释质量不均 | 对物理 TSS 去重并固定权重 | 不会自动修复错误注释或发现新 TSS |
| 生物重复少，批次差异明显 | 技术/生物重复/供体分层聚合，保留来源和覆盖 | 不能凭软件增加独立样本数 |
| 鸡等性染色体倍性不同 | 显式染色体倍性和 callability | 基因型、相位和 SV 输入质量 |
| 亲缘关系、群体结构使验证容易偏高 | 独立分组验证并防止重复区域泄漏 | 用户需正确声明家系/群体/实验分组 |
| RRBS 覆盖与 WGBS 不同 | 保留覆盖与缺失，不把未覆盖写成零甲基化 | 不同实验平台仍需匹配设计 |

适配来自这些明确的数据条件和检查，而不是给所有畜禽套用一组所谓最优权重。当前软件能提供可重复、可追溯的计算；跨组织、跨品种的生物学准确性需要在相应数据上验证。

## 13. 常见错误与继续阅读

- `PACE: command not found`：先 `conda activate pace`，确认安装和运行使用同一 Python 环境。
- 文件找不到：YAML 内路径相对于 YAML 所在目录；CLI 路径相对于当前目录；Docker 路径必须在挂载范围内。
- 输出目录已存在：使用新的输出目录，避免旧结果和新参数混在一起。
- 某基因全部 NA：检查固定活性组合、各正权重 TSS、零分母及 QC 原因，不要把 NA 改成 0。
- 数据合同不兼容：实际统一分辨率/尺度/窗口或重新拟合模型；修改字符串标识不能修复物理不一致。
- 自动参数仍为零：查看校准报告，可能是无标签、独立组不足、提升不稳定或本就没有独立增益。

[所有参数与默认值](parameters.md) · [字段字典](data_dictionary.md) ·
[输入准备](input_preparation.md) · [完整公式](FORMULA.zh-CN.md) ·
[训练与校准](training.md) · [局限与验证](limitations.md) ·
[故障排除](TROUBLESHOOTING.md)。
