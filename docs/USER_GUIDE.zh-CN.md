# PACE 完整使用说明

当前软件只支持以实测活性为基础的分析。本文按下载、输入准备、运行和结果解释展开。

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

需事先安装 Conda。该环境包含核心程序及常用基因组文件读取依赖。无需 GPU；[安装手册](INSTALLATION.md)给出了 Python venv、独立安装目录和 Docker 的安装方式。

Docker 从源码本地构建，容器入口已经是 PACE：

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  measured --config examples/measured/config.yaml --out results/docker_measured
```

### 1.3 运行合成示例

```bash
PACE demo --out results/demo
PACE validate --config examples/measured/config.yaml
PACE run --config examples/measured/config.yaml --out results/demo_measured
```

这些是合成软件检查，不能用作真实物种准确率。真实分析使用自己的实验数据并设 `execution_profile: research`。

## 2. 实测数据要求

活性必须来自合格 ATAC、DNase 或 H3K27ac。选择一个固定单层或双层组合。
接触默认采用实测；可显式配置适用的距离先验或接触收缩，结果必须保留来源。

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
| A_star | 固定检测组合的合格实测活性 |
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

[完整数学说明](FORMULA.zh-CN.md)给出总公式、完全展开式、接触策略、符号及零值/缺失规则。

基因表达权重没有再乘入主分数。同一基因的 TPM 权重同时乘到分子和分母会抵消；只在归一化后乘则改变分数含义。因此 RNA 作为独立注释或经过验证的机器学习特征输出。

## 4. 真实数据运行前需要哪些文件

软件从**已经处理的基因组数据**开始，不执行 FASTQ 比对、峰识别、变异检测、相位推断或通用坐标转换。所有文件必须使用同一物种、参考组装版本和染色体命名。

| 文件 | 是否需要 | 内容 |
|---|---|---|
| `units.tsv` | 必需 | 非重叠固定评分单元，默认 500 bp |
| `promoters.tsv` | 必需 | 基因、物理 TSS、方向、TSS 权重 |
| `candidates.tsv` | 必需 | 预先确定的元件—基因候选集合 |
| `samples.tsv` | 使用真实测量时必需 | 样本、个体、检测类型、生物/技术重复与背景 |
| `sources.tsv` | 有来源引用时必需 | 文件或 accession、处理方法、归一化、校验和 |
| `evidence.tsv` | 输入引用证据或导入解析结果时必需 | 证据身份和来源链；不能随意编造测量样本 |
| `observed_activity.tsv` | 必需 | 按固定单元汇总的定量活性及测量状态 |
| `observed_contacts.tsv` | 采用实测接触时 | 元件到各 TSS 的接触值、分辨率和尺度 |
| 接触先验目录 | 使用先验或收缩时 | 同尺度、同分辨率且适用背景的已拟合先验 |
| `expression.tsv`、甲基化、`features.tsv` | 可选 | RNA 和其他表观组学证据 |

各表必须是带表头的 UTF-8 TSV；不要把 Excel xlsx 直接重命名为 tsv。`NA` 表示缺失，0 表示真实测得的零。BED 和内部区间使用 0-based、左闭右开坐标；`tss0` 是 0-based 位点。[完整字段字典](data_dictionary.md)给出所有必需表头。

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

## 6. 实测数据运行配置与参数

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

## 7. RNA、甲基化及其他表观数据怎样使用

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

## 8. 如何看结果

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

## 9. 比较两只动物或两个处理

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

## 10. 为什么这个框架适配家养动物

| 家养动物分析中的实际问题 | PACE 的处理 | 仍需研究者解决的部分 |
|---|---|---|
| 物种、品种和参考版本多 | 显式绑定物种、组装、组织和模型范围 | 参考质量和跨版本坐标转换 |
| 表观检测不齐全 | 固定单层/双层活性组合，缺失保留 NA | 所选核心检测层仍需实测 |
| Hi-C 深度、分辨率和组织覆盖有限 | 接触测量合同检查、显式先验或可靠性收缩 | 适用先验拟合和独立评估 |
| 转录本/TSS 注释质量不均 | 对物理 TSS 去重并固定权重 | 不会自动修复错误注释或发现新 TSS |
| 生物重复少，批次差异明显 | 技术/生物重复/供体分层聚合，保留来源和覆盖 | 不能凭软件增加独立样本数 |
| 亲缘关系、群体结构使验证容易偏高 | 独立分组验证并防止重复区域泄漏 | 用户需正确声明家系/群体/实验分组 |
| RRBS 覆盖与 WGBS 不同 | 保留覆盖与缺失，不把未覆盖写成零甲基化 | 不同实验平台仍需匹配设计 |

适配来自这些明确的数据条件和检查，而不是给所有畜禽套用一组所谓最优权重。当前软件能提供可重复、可追溯的计算；跨组织、跨品种的生物学准确性需要在相应数据上验证。

## 11. 常见错误与继续阅读

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
