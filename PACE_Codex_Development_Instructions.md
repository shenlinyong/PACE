# 给 Codex 的开发指令：PACE 家养动物调控预测软件

修订日期：2026-09-17。将本文件完整交给Codex，并指定工作目录或现有仓库。本文件含完整数学附录，无需参考之前聊天。项目对外名称为 **PACE**，不加模型代际编号。机器所需的软件包版本、schema标识、注释版本和权重哈希仍保留。

## A. 任务、范围与交付原则

请直接实现可安装、可复现、可审查的Python科研软件，完成下列首次发布必交能力、测试、文档和GitHub仓库整理。不要只交框架或伪代码。先检查AGENTS.md、git status和已有项目，保护未提交修改、既有分析及API；沿用合理实现。本文是一套待实现的开发规范，不代表现有软件或性能。

三种模式必须完整保留：measured、hybrid、genome_only。主数学核相同，差异在证据解析。实现完成与真实权重训练、生物学验证分开验收。没有真实数据时，完成真实代码路径和明确标记的合成示例，不能训练随机示例权重后称支持任意家养动物。

### 首次发布必交

1. 统一数据校验、规范单元、主公式、多TSS、eta=0/1、接触观测/先验/一般收缩、来源与缺失输出。
2. 九组学标准表，以及BED/GTF、bigWig、cool/mcool、规范CpG计数、RNA表的适配。
3. atlas+启动子的WGS候选；定量序列模型适配、轻量参考CNN训练/推理、校准融合、SNV及可明确映射的小indel；输入中已报告而不支持的SV影响标记。
4. 三模式共同的bulk_proxy计算顺序；实际分母跟踪、共同分母比较、组成变化解释。
5. 有标签时的基础elastic-net分类器、分组调参、适当校准；无标签时独立运行公式模式。
6. 基础基准、重复稳定性、能力报告、离线demo、CI、文档和可上传仓库。

### 研究扩展，不阻塞首次发布

全基因组新候选扫描、任意组装/复杂CNV重构、copy_resolved_support、序列接触模型从零训练、Poisson–Gamma计数收缩、复杂后验/多层bootstrap、自动跨物种迁移。本文保留科学接口和边界，不要求首次软件全部解决。未实现扩展不得暴露返回成功的空CLI；请求时明确not_implemented。

不承担FASTQ比对、peak calling、WGS/SV calling、通用liftover或所有第三方caller格式。`.hic`不必原生支持，提供明确上游转换说明即可。前沿深度网络不要求重新实现。

当前任务是准备本地可审查仓库。只有当前执行会话已有明确远程仓库及发布授权时才推送；否则完成开发、检查和准确发布命令，不因缺远程地址停止实现，也不擅自公开数据。

## B. 软件结构与执行顺序

新项目可用发行名 `pace-livestock`、模块 `pace_livestock`、CLI `pace-livestock`；已有项目则保留其名称/API。名称可用性需实际检查。Python 3.11+，核心无GPU/网络依赖，IO、sequence、ml、dev分extras；采用实际测试的依赖范围，不虚构最佳版本。

```text
pyproject.toml
README.md / README.zh-CN.md
LICENSE / CHANGELOG.md / CONTRIBUTING.md
CITATION.cff 或有清楚待填字段的 CITATION.cff.example
.github/workflows/ci.yml
src/pace_livestock/
  cli.py / config.py / schemas.py / provenance.py
  catalog/         # 规范单元、候选边、标签映射
  io/              # tables, bed_gtf, bigwig, cooler, methylation, variants
  evidence/        # observed/predicted/prior/fused解析和适用范围
  core/            # activity, contact, tss, allocation, scoring
  sequence/        # 参考训练器、模型制品、个体序列准备
  learning/        # 分组、预处理、fit、calibrate、predict
  evaluation/      # compare、benchmark、replicate stability
  reporting.py
tests/
examples/measured/config.yaml
examples/hybrid/config.yaml
examples/genome_only/config.yaml
docs/model.md / data_dictionary.md / parameters.md
docs/input_preparation.md / validation.md / limitations.md
docs/github_release.md
```

实施顺序为：冻结schema与数学核→标准表三模式闭环→格式适配与真实序列计算→比较/学习/基准→安装、CI、文档。每阶段执行真实验收，不能到最后只剩目录。

必须有以下等价命令；名称如调整，同步所有帮助和示例：

```bash
pace-livestock validate --config examples/measured/config.yaml
pace-livestock capabilities --config examples/genome_only/config.yaml
pace-livestock prepare --config prepare.yaml --out prepared
pace-livestock run --config examples/measured/config.yaml --out results/measured
pace-livestock demo --regime measured --out results/demo_measured
pace-livestock demo --regime hybrid --out results/demo_hybrid
pace-livestock demo --regime genome_only --out results/demo_genome
pace-livestock fit-contact-prior --config contact_fit.yaml --out models/contact
pace-livestock train-sequence --config sequence_training.yaml --out models/sequence
pace-livestock predict-sequence --config sequence_prediction.yaml --out predictions
pace-livestock fit-fusion --config fusion.yaml --out models/fusion
pace-livestock prepare-genome --config individual.yaml --out prepared/individual
pace-livestock compare --config comparison.yaml --out results/comparison
pace-livestock variant-effects --config variants.yaml --out results/variants
pace-livestock stability --config stability.yaml --out results/stability
pace-livestock train --config learning.yaml --out models/classifier
pace-livestock predict-ml --config learning_prediction.yaml --out results/ml
pace-livestock benchmark --config benchmark.yaml --out results/benchmark
```

基础安装可运行measured；涉及序列网络才需要sequence extra。错误用非零退出码，已有目录不静默覆盖。示例运行/测试不下载数据或权重；首次安装依赖可能需要网络，不混淆两种“离线”。

## C. 输入与来源合同

UTF-8 TSV，统一NA表示，数值列类型严格。内部BED为0-based half-open，TSS为0-based碱基坐标；GTF转换保留版本。显式chrom alias表，组装不一致直接拒绝。配置未知键报错，路径相对配置文件解析。

### C1. 原始观测与模型估计必须分开

真实样本放samples.tsv；预测器和先验放evidence/models表，不能伪造一个“序列模型动物”。统一解析表使用 `evidence_id`，其中 `observation_sample_id` 允许空。技术或动物聚合用 `parent_evidence_ids` 及aggregation_method追踪组成，不当独立动物。

| 表 | 最少字段/主键 |
|---|---|
| `units.tsv` | element_id,chrom,start,end,anchor0,element_roles,canonical_catalog_id；ID/规范坐标唯一，默认完整非重叠500bp |
| `region_membership.tsv` | region_id,element_id,source_id,membership_rule；不让多来源重复增加单元 |
| `promoters.tsv` | gene_id,promoter_id,chrom,tss0,strand,pi,pi_source；gene+promoter唯一，一个物理promoter坐标一致 |
| `candidates.tsv` | element_id,gene_id,candidate_universe_id；边唯一，明确固定两个候选集合 |
| `samples.tsv` | sample_id,donor_id,assay,biological_replicate,technical_replicate,species,assembly,context_id,source_id |
| `observed_activity.tsv` | element_id,sample_id,assay,signal,measurement_status,callable_fraction,unit,normalization_id,window_id |
| `observed_contacts.tsv` | element_id,promoter_id,sample_id,contact_value,measurement_status,bin_pair_id,scale,resolution,source_id |
| `resolved_activity.tsv` | element_id,assay,observed_value,predicted_value,resolved_value,evidence_id,evidence_type,observation_sample_id,parent_evidence_ids,model_id,calibrator_id,fusion_weight,resolution_status,reason,unit,window_id |
| `resolved_contacts.tsv` | element_id,promoter_id,resolved_value,evidence_id,evidence_type,observation_sample_id,prior_id,reliability,resolved_mode,bin_pair_id,resolution_status,reason,scale |
| `features.tsv` | entity_type,entity_id,feature_name,value,evidence_id,status；E–G结构特征另含element_id,gene_id |
| `methylation.tsv` | chrom,dyad_start0,methylated_count,total_count,sample_id,assay；已合并正链CpG dyad，0<=m<=n |
| `expression.tsv` | gene_id,sample_id,tpm,status；transcript表另附到物理TSS映射 |
| `labels.tsv` | label_id,assayed_region_id,gene_id,context_id,perturbation_type,effect_direction,effect_size,label_status,assay_id,group_id,source_id；保留原实验单位 |
| `evidence.tsv` | evidence_id,evidence_type,source_id,parent_evidence_ids,model_id,unit,processing_method,checksum |
| `sources.tsv` | source_id,path_or_accession,source_type,assembly,processing_method,normalization_id,checksum |

`evidence_type=observed/sequence_prediction/contact_prior/fused/aggregate`。`measurement_status`描述实测：observed、unmeasured、low_coverage、unmappable、invalid、not_applicable。`resolution_status=resolved/unresolved/invalid`描述解析结果；物理删除另以structural_status记录，不能与NA混同。

所有观测行有真实sample；纯预测行sample为空但有evidence/model。缺失检测层无需制作全零文件。解析表由prepare/run生成或按相同schema导入，导入时仍完整校验来源/单位/模型范围。所有真正必需字段按regime和具体模式条件校验。

### C2. 规范单元与窗口

统一模式冻结非重叠参考网格，默认width=500、offset=0；atlas/峰/TSS映射到唯一单元。默认候选来源规则为至少1bp重叠，规则可配置并记录；跨来源去重。组装末端不完整单元标记并报告排除。TSS仍用精确坐标，不能取网格中心代替。

每个单元的activity读出、序列模型目标及监督标签窗口一致；输入window_id包含宽度、坐标规则和统计量。可视化用的合并区域不参与第二次评分。多个区域共用单元时，只在各自区域报告中复用，不能加回主分母。

旧区域表可用provided_regions兼容profile，但不同长度/重叠与信号统计明确保存；没有兼容预测器则只运行相应实测路径，不伪装成统一融合。转换表格本身不能凭空得到新网格定量值；需要重新从轨道提取或真实预计算表。

默认加入TSS所在单元并标promoter角色。候选半径示例5Mb，非物种最优常数；生成顺式稀疏邻域，不能构建稠密全基因组E×G。source candidate、actual normalization和comparison集合分别哈希。

### C3. 必需适配器

- BED/GTF：坐标检查、正负链TSS、物理去重、边界处理、规范单元及候选映射。
- bigWig：规范窗口内非负均值及callable_fraction；未存储区域只有来源明确表示已测零才按显式开关作0。负值轨道拒绝作为A。不混合均值/积分或不同规范化。
- cool/mcool：指定分辨率与raw/balanced、bin有效性、稀疏批量查询、共享bin_pair；元素锚点floor((start+end-1)/2)，TSS用tss0。近对角线显式选择先验/校正表/不可用。
- 甲基化：规范CpG计数、M_site/M_pooled及覆盖；未合并caller数据须有strand和参考才能转换。无参考CpG总数时比例NA。
- RNA：gene/transcript TPM及映射校验；gene TPM不平摊成所谓转录本实测。
- H3K4me1/H3K4me3/H3K27me3/CTCF：bigWig定量和可选BED overlap；CTCF锚点/区间内计数来自实测peak。无motif方向时方向NA，不从ChIP链推断。

## D. 三模式配置与资产依赖

以下为一个完整的synthetic measured配置合同；后两模式提供同一schema的覆盖项。开发时必须把它们展开为三个独立YAML并由CI实际读取，不能让用户手工猜合并方式。

```yaml
schema_version: 1
run_id: toy_measured
regime: measured
execution_profile: demonstration
estimand: bulk_proxy
target_level: individual
context:
  species: synthetic
  assembly: toy_assembly
  context_id: toy_tissue
inputs:
  units: units.tsv
  promoters: promoters.tsv
  candidates: candidates.tsv
  samples: samples.tsv
  observed_activity: observed_activity.tsv
  observed_contacts: observed_contacts.tsv
  evidence: evidence.tsv
  sources: sources.tsv
catalog:
  profile: canonical_grid
  width_bp: 500
  offset_bp: 0
  include_promoter_units: false  # 手算toy；生产参考配置为true
activity:
  panel: [ATAC, H3K27ac]
  combine: geometric_equal
  missing_policy: unresolved
contact:
  mode: observed
  scale: depth_normalized_contact
  prior_path: null
  near_diagonal_policy: prior_or_unresolved
  allow_prior_fallback: false
promoters:
  weights: provided
allocation:
  eta: 0
  missing_policy: fixed_gene_set
sequence:
  model_path: null
fusion:
  calibrator_path: null
genome:
  individual_id: null
  reference_path: null
  variant_path: null
  callability_path: null
  ploidy_path: null
  phase_policy: require_phase_or_single_variant_scenario
  unrecorded_site_policy: require_callable
multiomics:
  mode: annotate
  model_path: null
comparison:
  full_delta_requires_complete: true
  allow_conditional_intersection: true
  minimum_common_units: 2
output:
  format: tsv_gz
  retain_all_candidates: true
seed: 17
```

hybrid覆盖项如下；保留适用的真实观测输入，样本背景与calibration_target匹配：

```yaml
run_id: toy_hybrid
regime: hybrid
sequence:
  model_path: models/synthetic_sequence
fusion:
  calibrator_path: models/synthetic_fusion
genome:
  individual_id: synthetic_individual
  reference_path: genome.fa
  variant_path: sample.vcf.gz
  callability_path: callable.bed
  ploidy_path: ploidy.tsv
```

genome_only覆盖项如下；实验samples表为空，个体标识用genome.individual_id，需更多品种/家系信息时另设individuals表，不伪造实验观测：

```yaml
run_id: toy_genome
regime: genome_only
inputs:
  samples: null
  observed_activity: null
  observed_contacts: null
contact:
  mode: prior_only
  prior_path: models/synthetic_contact
  allow_prior_fallback: true
sequence:
  model_path: models/synthetic_sequence
fusion:
  calibrator_path: null
genome:
  individual_id: synthetic_individual
  reference_path: genome.fa
  variant_path: sample.vcf.gz
  callability_path: callable.bed
  ploidy_path: ploidy.tsv
```

外部模型已预计算结果可以导入，仍检查模型manifest与window/单位，不能以此绕过来源规则。没有WGS变异而仅提供参考序列时输出reference_context_prediction，不宣称个体化。

| 模式 | 条件必需 | 不必提供 |
|---|---|---|
| measured | 固定panel实测、候选/TSS、可用接触或先验 | 序列模型、VCF |
| hybrid | 可用实测及/或适用预测，明确来源；发生实际融合时须适用校准器 | 未测assay的空文件 |
| genome_only | 参考/个体序列、候选/TSS、适用定量模型、接触先验；个体化另需基因型可调用证据 | 新个体ATAC、Hi-C、RNA文件 |

每个run报告resolved_evidence_summary。hybrid退为纯实测或纯序列时必须说明实际使用路径，但不悄悄改panel或全局regime。

## E. 核心数值、接触和比较

严格按数学附录实现。函数不联网、不在predict时拟合。高动态范围用稳定log-sum-exp，真实零独立处理；不加通用epsilon。原尺度support溢出时可保留log_support与有效分数，原尺度NA并标overflow。

对已知完整联系集合全零的元件，eta=1虽B无定义，其全部支持明确为0并保留零状态；未知联系不可当此情况。eta=0不要求B。全基因分母零/空返回NA和原因。

fit-contact-prior在指定同物种训练区域用有效bin pair距离分箱，均值包含零。首次实现对正均值箱做明确的log线性加权拟合，保存每箱有效pair数、zero-bin处理、参数范围和留出残差；权重不是最优统计声明。a>0、gamma>0，否则报告拟合失败或不符合该先验族。near-diagonal使用d_min平台是独立规则，不能把它藏在拟合外推中。默认不对跨染色体联系套幂律。

一般收缩输入Cobs、Cprior、r和r来源，单位一致；缺少校准参数时允许observed，不假称做了收缩。计数Poisson–Gamma留作研究扩展，不阻塞核心。

compare不能inner-join原分数后直接相减。先核对estimand、target_level、catalog、panel、pi规则、单位、B候选集合及参数；再从support按共同可评分单元重算各侧分母。至少2个共同单元且双方分母正才输出conditional组成差异；完整Delta还要求计划宇宙无必要技术NA和结构关系可解释。不要修改B的候选基因集合来改善覆盖。

输出original_score、common_score、comparison_universe_id、full_delta/conditional_delta和reason；原始全覆盖结果另报。Delta support只有共同尺度成立才计算。差异分析同时报告Delta A、gene total_support和分母变化，不从Delta PACE推表达效应方向。

stability仅做样本/动物层面的一致性，技术重复不当独立动物，常数向量相关为NA。缺失不同按共同分母重算并标条件。首次实现不要求小样本输出置信区间；有研究重采样模块时需按独立单位和配对结构抽样。

## F. 定量序列、WGS与融合

### F1. 最小可训练序列模型

sequence extra提供实际train/save/load/predict。默认输入8192bp one-hot DNA，中央500bp规范单元为标签目标，输出species/context/assay头的非负信号。N作为未知掩码，N比例超过配置阈值则不可用。

可实现参考CNN：Conv1d通道64/128/128、kernel15/7/5、GELU、pool4/4/4；最后按特征bin与中央目标区的重叠作局部池化，可拼接全窗口上下文，再线性头+Softplus。保存准确padding/坐标映射；不能只做全局池化丢失目标位置。架构是基线，可经独立验证调整。

masked Huber作用于log1p(pred/s)和log1p(label/s)，每个head的s从训练正信号确定并冻结，全零/无标签头不可用。反向互补增强用于这些非链特异活动任务。默认AdamW、lr=1e-3、batch64、50轮上限、patience5、seed17均为工程起点。模型不自动获得精确后验或每拷贝剂量单位。

训练标签和输入窗口按染色体/区块隔离，重复、同一区域其他组织及增强副本不能跨fold泄漏。外部模型manifest至少含模型ID、hash、species/assembly/context、assay、input_length、output_window、output_type、signal_unit、normalization_id、target_level、training/calibration/test来源、真实验证报告路径与is_synthetic。

peak_probability/logit不能直接进入定量融合，分类概率校准不解决量纲。normalized_signal仍要检查单位与窗口。外部contact_enrichment须配距离期望还原，首次发布可完全使用距离先验。

### F2. WGS范围与可调用性

支持规范SNV和明确可映射的小indel；REF校验、多等位、重叠冲突、GT/phase_set、callability和ploidy有独立测试。未记录VCF位点默认只有在callable区域才按参考处理；没有callability时仅显式研究假设允许，输出unknown/assumed_reference比例。缺失GT不可变成0/0。

相位不确定时默认不给虚构完整单倍型，提供单变异情景或明确phase scenarios。跨phase block没有相位关系时也要说明。reference-only输出不得命名为个体基因组结果。

固定bulk_proxy主路径先逐assay平均单倍型定量预测，再融合和构造几何活性；接触也先解析到同一边际层级。固定二倍体rho=1/2只是声明的可加平均proxy；不能把CNV剂量和bulk测量直接套进去。copy_resolved_support不是首发主输出。

indel后的目标窗口必须可说明地对应模型训练目标；目标区长度、边界/TSS或对应关系不再匹配时，该完整PACE比较不可用。局部序列差可以单列，但不能硬填到主分数。明确记录哪些小indel路径已支持，哪些因窗口合同返回unsupported，不能笼统声称支持全部indel。

仅识别用户VCF/映射中已报告的SV/CNV，记录structural_variant_source。标记受影响模型输入窗口、E/TSS和涉及距离关系。没有SV来源时写not_assessed，不声称不存在SV。无法解析的影响边不继续用参考坐标评分。复杂拷贝重构、任意组装映射是研究扩展。

### F3. 校准器

fit-fusion按附录凸最小二乘实现，保留质量层、尺度、样本量、来源、目标估计层级和独立评价。individual_state必须针对同动物适当独立的测量，不借其他动物的差异当技术噪声；population_mean的校准器不能用于“恢复个体状态”。无适用校准时单来源选择，不猜融合比例。

w=0/1分支不要求未使用来源；D=sum[(zobs-zseq)^2]=0时不可识别。最低必要样本量/质量层合并规则由训练协议明确，不能只用一条边就称已校准。参考训练器和融合器均有真实代码，CI不要求训练出高生物学准确率。

## G. 九组学学习与标签

九组学全部具备实际提取、状态和导出接口。默认附加标记annotation_only，适用已训练分类器才active_ml。它们不全乘进主公式。

特征schema冻结：log1p(A_used/Cbar/distance)、主分数、可用E/P标记、CpG及覆盖、CTCF结构、缺失指示、regime与来源。gene TPM默认关闭，仅作独立消融。promoter特征按pi汇总；正pi启动子缺失时该特征NA及missing_pi_mass，不重新分配pi。

固定变换→预定义少量交互→训练折median填补和median/IQR缩放；IQR0用尺度1并记录，整列无训练信息则剔除。核心A/C/score无效不能靠填补救回。分类器范围须覆盖实际regime和预测来源；超范围的score/概率字段不能冒充验证。

实现elastic-net logistic及base-only对照，截距不罚，样本权重默认1。若使用库参数C/l1_ratio，文档写清和目标损失的比例关系。内部按可行分组选择超参数，外部固定test保留；独立校准集或严格OOF用于sigmoid calibration。没有足够分组/类别就明确不能调参/评价，可接受预先固定参数研究运行，不造样本补折。

阳性方向、效应量和功效写入labels规范。默认正增强联系为抑制元件后基因下降；上升单列potential_repressive_or_complex。未测试或低功效不作负例。区域标签不得复制到多个tile冒充多个独立阳性。首次训练仅允许明确一对一映射，模糊映射另列；区域级评价按冻结的单元组聚合支持并记录分辨率。

预测模型保存线性系数、预处理、缺失列、类别和校准参数，优先JSON/NPZ；序列权重使用安全张量制品配JSON。不得在加载模型时重新拟合。无标签时公式路径不受阻；无校准概率字段NA，不平均ML与PACE。

## H. 基准、能力报告与输出

capabilities报告：implementation_status、weights_status(absent/synthetic/real)、task_validation_status、scope、execution_profile和每个能力的阻塞原因。不能只设置一个validated=true开关替代真实报告。

demonstration仅用合成数据/权重；research允许真实但未完成特定任务验证的模型并明确状态；validated核对相应task/context的证据。research/validated不能使用synthetic权重，demonstration可以通过run或demo调用。WGS权重缺失不代表软件实现失败，必须准确说明缺的是模型资产。

run至少输出：

| 文件 | 内容 |
|---|---|
| scores.tsv.gz | 全候选，包括NA及原因、A/C/B/support/denominator/score、normalization_status、来源、估计对象 |
| gene_summary.tsv | 候选数、可评分数、正支持数、条目覆盖、分母、归一化集合ID |
| resolved_activity.tsv / resolved_contacts.tsv | 原观测、预测、融合、来源、模型、窗口、单位和解析状态 |
| multiomics_features.tsv.gz | 九组学特征、覆盖、角色和来源 |
| qc_report.json / report.md | 数据问题、模型范围、能力状态、未支持变异、解释限制 |
| resolved_config.yaml / run_manifest.json | 实际参数、schema/软件/权重/输入哈希、seed、环境、各宇宙ID |
| comparison.tsv（如请求） | 原分数、共同分母分数、完整/条件Delta、A/S/基因总S、比较原因和ID |

完整评分集合合法时每基因分数和为1。partial也可能和为1，因此必须带normalization_status；单个正支持得分1不算强证据。B在eta0时可NA。未请求学习、不确定性或未有资产时附加字段NA。

benchmark读取真实可用分数或生成明确配置的基础对照。首发基础对照包括负距离排名、ABC-style单TSS、PACE eta0/eta1；单TSS选择规则预先冻结，不能用测试表现选启动子。gABC/外部方法可接受带版本和配置的分数表，不假装本软件已完整复现外部包。缺标签/权重的方法标not_available，不填0，也不阻塞其他合法对照。

默认报告average_precision(AP，说明定义，不混同梯形PR面积)、可评估的AUROC、冻结阈值precision/recall、覆盖和各失败原因。测试PR曲线上的固定召回precision可描述，但不作为部署阈值来源。区域/基因/距离/表达/物种/检测层分层样本数同时报告。

共同候选评分比较要重算共同分母；端到端另报候选发现及不可评分漏检，在预定完整测试范围的阳性上计算召回。不要把缺预测值填0后当作原始连续得分；用coverage和missed-positive计数计算端到端覆盖/召回，候选内曲线单列。

序列位点预测、个体Delta、功能联系三个任务分别验证，QTL关联仅作辅助支持。没有真实标签只做软件验证报告，不能生成“提高多少”的生物学结果。

## I. 必须通过的验收

数学预期详见附录，测试使用独立期望常量，不把实现复制一遍当断言。

1. 核心三元件例、活性4×9、多TSS、eta跳过B、真零/NA、全零分母、共享bin、启动子去重均正确。
2. bulk反例必须输出5/9；研究copy求和反例3/7不能混入主模式。三个regime使用相同aggregation代码。
3. 分母不同不能直接Delta；共同集合须从support重新归一化，B集合保持原定；完整Delta与conditional分开。
4. 规范单元去重；同一峰重复来源不增加支持；任意输入顺序不改变catalog和结果；输出500bp模型拒绝不匹配长区域。
5. 观测/预测证据样本字段正确：genome-only无假样本、无空实测文件依赖；纯先验来源可追踪。
6. 融合sqrt(40)-1，w=0/1和D=0边界；individual/population校准器错用拒绝；概率与强度类型错误拒绝。
7. VCF REF错误、GT缺失、未记录不可调用位置、相位区块、倍性、indel窗口不匹配和已报告不支持SV分别得到明确状态。
8. CpG (1,2),(8,8)得到0.75/0.9；双链不重复、无CpG与零甲基化不同。
9. 标签上升/低功效/未测试不进入普通二分类阴性；多tile标签不复制；测试集极端值不改变训练预处理。
10. CNN训练冒烟测试仅验证loss/梯度有限、参数确实更新、mask正确、保存加载一致；确定性固定权重fixture检查序列变动/中央定位/三模式闭环。不设置“几个epoch必须达到某准确率”的脆弱CI门槛。
11. 注释层变化不影响主分数；ML存在时确实读取其特征；缺额外组学不破坏核心。
12. 从wheel干净安装并跑三模式demo；基础/IO/sequence/ML extras分job，运行不下载数据权重；execution_profile!=demonstration时拒绝合成权重；run读取demonstration配置时正常允许，不能按CLI命令名称拦截。
13. 测试候选/矩阵算法不构建全基因组稠密E×G。记录中等合成规模实测耗时/内存及环境，不承诺未测性能。

数据字典与示例schema必须一致。float64手算可用atol=1e-10；不要对无关文档格式写无意义测试。实际未执行的测试明确列出原因，不声称通过。

## J. GitHub与最终交付

中英文README讲清最少输入、三个模式、公式含义、真实模型资产需求和不支持范围。不要将“支持WGS输入”写成“已验证所有家养动物个体预测”。

新代码可用MIT，现有仓库保留许可证，引用第三方实现前核查条件。作者只用用户已确认的身份和项目资料；没有完整作者表时不编造其他作者/DOI。CITATION可用明确模板，核心实现不因元数据不全停止。

GitHub只带代码、文档和小型synthetic fixtures；真实测序数据、私有路径、凭据、大权重/结果不自动提交。CI在Linux完成基础、extras、wheel与demo，不依赖私有账号。可扩展其他平台，但不阻塞首发。

对外标题、README、图表统一称PACE。内部包管理需要合法初始版本，schema/模型文件和注释需要追踪标识；这些不拼进模型名称。不重写历史、不force push、不擅自正式release。若已有明确发布授权，按授权完成；否则交付准确远程配置步骤。

最后报告实际实现能力、文件、安装/运行命令、真实测试结果、权重可用性、验证状态、局限及仓库状态。不得只重复计划。

## K. 完整数学与科研合同

以下是本开发指令的组成部分。整理进docs/model.md、parameters.md、validation.md和limitations.md，并保持公式、字段、默认值和边界行为一致。

<!-- CANONICAL_MODEL_BEGIN -->

# PACE：家养动物增强子—基因调控预测模型

修订日期：2026-09-17。对外名称统一为 **PACE**，不附加模型代际编号。本文定义模型、数据条件和验证要求；尚不代表软件已经实现、目标物种权重已经训练或性能提升已经成立。配置格式、软件包、注释和权重文件仍须保留内部追踪标识。

## 1. 要解决的问题与结果含义

PACE 估计某个候选调控单元 E 对基因 G 的相对增强支持。主分数由活性、接触及可选的跨基因分配构成。它不是基因表达贡献百分比、因果概率或表达变化倍数，也不直接预测抑制性调控强度。

主输出统一使用 `estimand=bulk_proxy`：表示一个明确样本或群体背景下，由各检测层的边际信号构成的调控支持代理。`target_level=individual` 与 `population_mean` 另行区分，不把两者的训练误差、融合系数或验证结论混用。bulk_proxy 不是单倍型物理支持总和；它与普通 bulk 组学的测量层级一致。

| 数据模式 | 活性来源 | 接触来源 | 输出定位 |
|---|---|---|---|
| measured | 合格实测 | 实测；或显式同物种先验 | 实测支持的相对调控预测 |
| hybrid | 合格实测、适用定量序列预测及经过校准的融合 | 合格实测和可用先验 | 混合证据预测 |
| genome_only | 适用物种、组织的已训练定量序列模型 | 同物种距离先验；经过验证的序列接触模型可作扩展 | 指定组织下的遗传调控潜能 |

新个体无表观数据时，可以使用同物种既有模型；整个目标物种缺少功能数据时，只能显式进行跨物种探索。软件接口支持某物种，不等于已有该物种的有效权重。DNA 无法唯一恢复个体当时的营养、感染、激素、细胞比例和表观记忆。

`pace_score`、独立学习分数 `pace_ml_score`、经适用校准的 `pace_ml_probability` 分列；不将公式分数和分类概率平均。证据来源、缺失、覆盖、适用范围另外输出。

## 2. 总公式与计算顺序

省略物种、组织/状态、个体或群体背景下标。一次运行必须明确这些条件。

\[
S(E,G)=A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta},
\qquad
\boxed{\mathrm{PACE}(E,G)=\frac{S(E,G)}{\displaystyle\sum_{e\in\mathcal E^{score}(G)}S(e,G)}}.
\]

\[
A_\star(E)=\prod_{m\in\mathcal M}x_{\star,m}(E)^{1/|\mathcal M|},\qquad
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),
\]
\[
\widetilde C(E,t)=r(E,t)C_{obs}(E,t)+[1-r(E,t)]C_{prior}(E,t),\qquad
B(E,G)=\frac{\overline C(E,G)}{\displaystyle\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
\]

| 符号 | 定义 |
|---|---|
| E | 规范评分单元，可对应一个增强子或增强子的一部分；展示区域与评分单元分开 |
| G、t | 目标基因、去重后的可信物理启动子/TSS |
| M | 本次固定活性检测层组合：ATAC、DNase、H3K27ac、ATAC+H3K27ac、DNase+H3K27ac 五选一 |
| x_star,m | 与模型单位、窗口、背景匹配的非负解析信号，来源可为实测、序列或融合 |
| pi | 同一基因启动子权重，非负、和为1，不依赖E |
| r | 接触观测权重，范围[0,1]，来源明确 |
| B | 同一元件在固定候选基因集合中的接触分配，不是边界强度或真实资源守恒 |
| eta | 首次实现限定0或1；默认0，1用于预先指定的对照 |
| E_score(G) | 原定候选宇宙中，本模式能够计算支持的单元集合 |
| G(E) | 固定候选基因集合，不随本次分数或表达阈值临时改变 |

三个模式的主流程相同：**统一评分单元与测量层级 → 逐检测层解析/融合 → 活性几何均值 → 多TSS接触 → 可选B → 每基因归一化。**不能在某模式改成先计算单倍型支持再相加。

eta=0直接跳过B，避免0^0。r=0/1时只要求实际使用的来源，不能让0×NA污染计算。所有必要信号合法且所有候选接触均为0时，B可为NA，但该元件支持明确为0；不能误报为技术缺失。接触或必要输入未知则不同，按不可评分处理。每基因支持总和为0时分数NA，不返回均匀分数，也不加任意epsilon或未知支持U_G。

## 3. 规范评分单元：避免重叠和长度造成重复计数

统一模式默认使用冻结、非重叠的参考坐标网格，启动值 `unit_width_bp=500`、`grid_offset_bp=0`。单元以0-based half-open坐标定义。atlas、实测峰、启动子和未来扫描结果只决定哪些单元被纳入及其来源，不把同一单元重复加入分母。组装末端不足一个完整单元时标记并排除，记录覆盖损失。

500 bp是工程起点，不是家养动物的最优生物学常数。改变宽度/网格偏移时重建catalog、重做匹配信号与序列标签，并作为独立敏感性分析。默认不把多个重叠500 bp预测窗口合并成更长区域后，仍用一个窗口的预测代表整个区域。

每个规范单元的活动统计、序列目标窗口和标签必须一致。峰区域可以合并用于展示；区域级相对支持只能对其所含、不重叠单元的支持求和，不再次把“区域”作为另一个元件加入同一分母。区域分组规则须冻结，不能据评分结果挑有利的组合。

保留 `provided_regions` 兼容模式以接收已有区域表，但其长度/重叠处理必须明示，输出不同的 `scoring_profile_id`；缺少窗口匹配的预测器时不允许在该模式做定量序列融合。不宣称不同评分单元方案的分数可直接比较。

默认加入含可信TSS的规范单元，并标记promoter角色，不重复添加与已有单元重叠的启动子元素。不强制其分数为1。仅远端模式单独标记，不复用主配置阈值。TSS本身仍是精确坐标，不被网格起点替代。

`canonical_catalog_id`标记全部规范单元及定义；`candidate_universe_id`标记冻结的E–G边；`normalization_universe_id`标记每个基因本次实际分母中的单元。三者不同。

## 4. 活性、重复与测量层级

双层活性为sqrt(可及性×H3K27ac)，单层直接使用。ATAC和DNase不能作为两份独立可及性层同时相乘。固定panel后，不因某个元件缺一层而偷偷改成单层；可以另建明确的单层运行。

输入须是单位、规范化协议、统计量和窗口明确的非负信号。负的log fold-change不能直接进入几何均值；不静默裁剪。真零保留，未测/低覆盖/不可比对为NA。技术重复原始计数先合并再规范化；已规范化轨道使用明确的加权均值，不直接相加。动物层面重复等权汇总时，结果是群体平均代理，不是每只动物结果。

三模式的主估计都先得到相同层级的各层信号。对于固定倍性、无影响位点对应关系的结构变异的序列输入，可显式采用对称可加的平均代理：
\[
x_{seq,m}^{bulk}(E)=\sum_h\rho_h x_{seq,m,h}(E),\qquad
C_{prior}^{bulk}(E,t)=\sum_h\rho_h C_{prior,h}(E,t),\quad\sum_h\rho_h=1.
\]
固定二倍体可用rho=(1/2,1/2)，单倍体用1；这表示所声明的平均代理，不是证明染色质信号按剂量线性变化。rho不能用测序深度或未验证的等位表达比例替代。单倍型局部联系不同时，这种边际代理仍不等于平均物理支持，必须保留解释边界。

bulk实测不能复制给两个单倍型，也不能默认除以2称为等位实测。只有可靠等位信号、接触和每拷贝单位足够时，才可另做 `estimand=copy_resolved_support` 研究分析：先算各拷贝支持、相加再归一化。它与主bulk_proxy是不同对象，输出单独命名；首次发布不以该扩展为必需能力。

## 5. 活性融合与校准

将单位匹配的观测和序列预测变换到相同空间：
\[
z_{obs,m}=\log(1+x_{obs,m}/s_m),\qquad z_{seq,m}=\log(1+x_{seq,m}/s_m),
\]
\[
z_{\star,m}=w_m z_{obs,m}+(1-w_m)z_{seq,m},\qquad
x_{\star,m}=s_m\operatorname{expm1}(z_{\star,m}).
\]

s_m>0来自训练/校准协议并冻结。w_m是观测融合权重，不是B，也不是活性层之间的权重。measured取w=1；genome_only取w=0；hybrid有适用校准器才融合。无校准器时优先合格实测，缺观测且有适用预测时用单来源预测，保留实际路径；两者无效则NA。不设通用0.5。

严重污染、错组织或明显不相容批次不得只给小权重混入。序列输出若是峰概率/logit，即使已校准分类概率，也不是定量强度；只有额外的定量校准通过检验后才能进入x_seq。

按预定义质量层b拟合凸组合：
\[
\widehat w_b=\arg\min_{0\le w\le1}\sum_{i\in b}[z_i^*-(wz_{obs,i}+(1-w)z_{seq,i})]^2.
\]
闭式解为clip(sum[(z_obs-z_seq)(z*-z_seq)]/sum[(z_obs-z_seq)^2],0,1)。分母0表示不可识别，不能制造精度；按明确单来源规则输出。

校准器必须声明 `calibration_target=individual_state` 或 `population_mean`。前者需要同动物的适当独立测量或有说明的独立技术拆分；用另一只动物作真值会把真实个体差异当噪声。后者可用动物间共识，但只校准群体平均。普通序列模型不能因此恢复个体非遗传状态。

序列训练、s_m、质量分层、候选选择阈值、融合权重及最终测试的角色全部记录。低深度输入不能同时混入所谓独立高深度真值；优先独立重复，拆分读段则记录依赖及局限。校准评价至少包括低深度实测、序列单来源、融合和独立参考，并同时报告覆盖。

独立、近似无偏且误差已估计时，w=v_seq/(v_obs+v_seq)可作为研究估计；相关误差还需协方差。首次发布默认使用上面的直接校准，避免把模型集成分歧或FRiP当成完整测量误差。融合是点估计，不自动构成后验。

## 6. 接触、多TSS与分配

接触分为observed、prior_only、shrinkage。输入必须是非负、包含距离背景且尺度可比的量；O/E或log O/E须有匹配期望才能还原。环调用P值、相关系数不可直接当C。

cooler稀疏矩阵的未存储pixel，仅在有效bin和明确存储语义下才代表0；无效bin不是0。同一bin pair的多条E–TSS边共享测量和抽样标识。实测联系、参考组织联系与先验分别标记；CTCF不能代替Hi-C或在主公式外再奖励一次。

可用同物种距离先验：
\[
f_s(d)=a_s\left(\frac{\max(d,d_{min})}{d_{ref}}\right)^{-\gamma_s},\quad a_s>0,\ \gamma_s>0.
\]
拟合使用有效bin pair，包括真实零；保存数据尺度、分辨率、距离范围、训练区域、组织和留出残差。d_min是近距离规则，不能与“加一个5000 bp伪计数”混淆。首次实现采用清楚的分箱拟合，不声称其误差权重是最优；不能仅拟合正pixel。近对角线策略单独记录，默认使用同尺度先验的近距离平台；这是建模选择，需敏感性检查。没有先验或有出处的校正时不可评分并报告启动子分母缺失。

普通收缩使用有出处的r。原始计数Poisson–Gamma为研究扩展：
\[
N\mid\lambda\sim\mathrm{Poisson}(L\lambda),\quad
\lambda\sim\mathrm{Gamma}(\kappa\mu,\kappa),\quad
\widetilde C=(N+\kappa\mu)/(L+\kappa).
\]
Gamma第二参数为rate。L、kappa、mu均正；L不能取本边count，必须处理所声明的曝光/偏差模型。ICE/KR平衡值不是Poisson整数计数。不支持该模型时不得假装已进行贝叶斯校正。

TSS以gene、chrom、tss0、strand精确去重；多个转录本同TSS只算一次。GTF正链tss0=start-1，负链tss0=end-1。不同基因共享启动子标歧义。pi优先匹配CAGE/RAMPAGE；无可靠起始信息时对可信TSS等权。gene TPM不提供TSS使用率，H3K4me3也不直接当起始次数。pi>0的必要接触缺失时不删除该TSS重分配pi；pi=0不参与求和。

E–G distance_bp定义为元件锚点到可信TSS的最小绝对距离，接触仍按每个E–TSS距离计算。默认仅顺式联系。

B分母始终使用原定G(E)。eta=1若有必要联系未知，B不可评分，影响该E对全部相关基因的结果；完整报告这类覆盖损失。不能因某个基因缺数据而缩小B分母。B可能惩罚真实共享增强子，因此默认eta=0。全部已知接触均0时支持0的规则见第2节。

## 7. 缺失、实际分母与跨样本比较

至少区分三种状态：

| 情况 | 数值处理 |
|---|---|
| 确证物理删除 | 在拷贝身份和目标基因仍可解释时，相关物理单元支持为结构性0；主bulk涉及未校准剂量或结构关系时标不支持，不擅自推总量 |
| 技术不可测、映射失败、必要证据缺失 | NA，不能当遗传效应 |
| 未过峰/候选发现阈值 | DNA区域仍存在；比较时在共同坐标重新定量，不能置0 |

`normalization_status=complete/partial/zero_support/empty`与分数一起输出。complete只表示预定义候选均被处理，不表示发现了所有生物学增强子。partial分数条件于可评分子集。条目覆盖率不能解释为捕获的真实调控质量。

比较样本、模式或eta时，相同候选文件不保证分母相同。必须从各方法的未归一化支持出发，在共同可评分集合重新归一化：
\[
\mathcal I_G=\bigcap_k\mathcal E_k^{score}(G),\qquad
P_k^{common}(E,G)=S_k(E,G)/\sum_{e\in\mathcal I_G}S_k(e,G).
\]
B仍按各运行同一原定G(E)计算，不为了比较删基因。保存comparison_universe_id，原分数不覆盖。交集不足2个单元、某侧分母0或无法建立同源对应时，组成差异不具有所需比较信息，返回NA及原因。

完整个体Delta要求同一计划宇宙、单位、估计对象、参数和可解释的结构状态；任何必要技术NA使完整Delta为NA。可以另报共同可测集合上的 `conditional_delta_pace`，明确其条件，不能称完整个体遗传效应。新增/删除结构与仅越过活性阈值分别报告。

同时输出Delta A、Delta support、基因总支持与Delta PACE；没有共同尺度时绝对支持差NA。相对份额下降不一定是该元件活性下降。所有候选同时增大，PACE也可能不变。

## 8. 九类组学与独立学习扩展

| 数据 | 必须保留的量 | 默认用途 | 有适用功能标签时 |
|---|---|---|---|
| RNA-seq | gene/transcript TPM、检测状态 | 表达注释，不乘主分数 | 可选表达特征，另做表达分层 |
| ATAC/DNase | 规范单元信号、可测比例 | 活性 | 基础特征及预定义交互 |
| H3K27ac | 信号、单位、input处理信息 | 活性 | 基础特征及状态交互 |
| H3K4me1 | E信号、peak overlap | 元件状态注释 | 特征 |
| H3K4me3 | TSS窗口信号 | 启动子状态注释 | 特征 |
| H3K27me3 | E/P信号、broad-peak overlap | 状态注释 | 特征，不强制固定负系数 |
| CTCF ChIP-seq | E/P占据、区间位点数/强度 | 结构注释 | 有motif才用方向特征 |
| DNA methylation | E/P甲基化、CpG及覆盖 | 状态注释 | 分区域且带覆盖的特征 |
| Hi-C | 接触、尺度、分辨率、bin状态 | 接触 | 基础特征 |

每层标明active_core、annotation_only、active_ml或disabled，并另标observed/predicted/fused。缺少的额外组学不生成假数据。H3K4me1不能简单等同活性输出，甲基化也并非在所有位置统一抑制转录因子结合。[Dorighi等](https://doi.org/10.1016/j.molcel.2017.04.018)、[Hu等](https://elifesciences.org/articles/00726)

CpG计数m_j、n_j先统一到正链dyad坐标，禁止重复合并。定义M_pooled=sum(m)/sum(n)，M_site=mean(m_j/n_j)。后者默认作为特征，前者也输出；记录覆盖和CpG数。无CpG、未覆盖、低覆盖和真实零不同。无参考CpG总数时覆盖比例NA；常规bisulfite通常不区分5mC/5hmC。

学习模型使用固定特征映射的elastic-net logistic：
\[
p=\sigma(\beta_0+\boldsymbol\beta^T\phi),\qquad
\mathcal L=-\frac{\sum_i\omega_i[y_i\log p_i+(1-y_i)\log(1-p_i)]}{\sum_i\omega_i}
+\lambda_1\|\beta\|_1+\frac{\lambda_2}{2}\|\beta\|_2^2.
\]
omega默认1，非负、总和正；lambda非负，截距不惩罚。phi包括log1p(A)、log1p(C)、log1p(distance)、PACE和可用额外特征。用omega区别于活性融合w。仅核心可评分边可训练/推理，核心缺失不得靠中位数填补。额外特征可用训练折中位数和缺失指示；median/IQR缩放及所有交互定义在训练侧冻结。

默认阳性是抑制增强子后靶基因下降，满足预设效应量、显著性和功效条件。表达上升标为potential_repressive_or_complex；未测试、功效不足、明显存活混杂或基因无法检测的不作普通阴性。区域扰动覆盖多个规范单元时不能给每个单元复制同一阳性标签；首次实现只用可靠一对一映射训练，模糊映射另报。区域级评价可按预注册单元分组合并支持，不宣称解析了区域内每个单元的因果作用。

预测特征来源与regime须在分类器训练/校准范围内。无功能标签时公式仍运行；无独立概率校准时只报classifier score。经校准概率仅适用于记录的实验及候选抽样分布，不是普适因果概率。需要base-only分类器作为对照，检验额外组学的增益。分类系数不等于各标记的生物学因果效应。

基因总表达量R(G)>0在分子分母统一相乘会约掉，R=0则0/0。归一化后再乘会改变分数定义；因此RNA默认注释而非主公式乘数。

## 9. WGS、序列模型与模型资产

首次发布必须支持同物种atlas+可信启动子的genome_only模式，不要求新个体ATAC峰。这能给已有候选评分，不能声称发现了新品种全部增强子。全基因组扫描及其发现阈值属于研究扩展，另评估端到端漏检。

输入为参考FASTA、注释、VCF/BCF或规范化个体序列、callability、倍性及明确组织。VCF无记录不自动意味着可靠参考基因型；优先gVCF或callable BED。只有显式研究假设才能把未记录位置当参考，并保存不确定范围。显式缺失GT不得当0/0。

首次实现SNV及规范化小indel，默认已测试长度不超过50 bp；校验REF、多等位、重叠冲突、相位区块和模型输入窗口。未相位杂合不能拼任意单倍型；可做有明确假设的单变异REF/ALT情景，不能称已重建个体完整基因组。

indel会改变映射和距离。主bulk比较只接受评分目标、TSS和窗口定义仍能对应的情形；不能将参考500 bp区间变为不同长度后仍冒充相同标签窗口。参考与个体输入要保留output_target_id和映射状态；边界/元件/TSS被破坏或目标定义不能一致时，该比较NA。模型不支持的情况下，局部variant情景可另输出序列信号差，但不强行升级为完整PACE差值。

仅识别输入中已报告的SV/拷贝变化，标记受影响序列窗口、E/TSS及距离关系；本软件不从普通SNV VCF发现遗漏的SV，也不承担WGS calling。复杂重排、CNV剂量重构和任意组装映射不作为首发必交。鸡性染色体等倍性须显式输入，不全设二倍体。

提供一个轻量定量CNN参考训练器及外部模型适配协议。启动输入8192 bp、中央500 bp输出目标，保留中央位置特征，可拼接全窗口上下文；训练标签是与推理一致的规范化信号。它是可复现基线，不是最佳架构承诺。masked Huber(log1p(pred/s),log1p(label/s))可用作训练目标；s只在训练数据确定。

输出模型卡分开记录软件能力、权重来源、适用范围、跨位点验证、个体差异验证和联系验证；布尔true必须指向实际验证报告。DeepFARM是家养动物序列分类的相关先例，其峰分类输出不能直接当这里的定量活性。[DeepFARM论文](https://doi.org/10.1093/nargab/lqaf139)

| 执行环境 | 允许权重 | 解释 |
|---|---|---|
| demonstration | 明确标记的合成/固定权重 | 测试代码，不证明生物学能力 |
| research | 真实训练、适用或明确外推的权重 | 可计算但可能未完成目标任务验证 |
| validated | 真实权重且指定任务/背景有验证证据 | 只在相应范围使用验证表述 |

无真实权重时，三模式软件仍可完成实现与合成测试，但不能宣称安装后即可可靠预测任意家养动物。WGS不能替代训练阶段所需的功能数据。

## 10. 科研假设与最小验证矩阵

主要待检验假设：①缺失/降深时校准融合优于最佳单来源；②适用物种的接触和TSS处理在相同候选下改善性能或覆盖；③个体序列在独立动物或等位数据上提供额外信息。多TSS、ABC乘积和基因分配已有相关工作，不能单独作为原创主张。[ABC官方方法](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/methods.html)、[gABC/STARE](https://github.com/SchulzLab/STARE)

| 要验证的能力 | 必需证据 | 对照与报告 |
|---|---|---|
| 数学/软件正确 | 手算和合成输入 | 边界状态、可复现性；不报生物学准确率 |
| 跨位点活性预测 | 留出区域及适用动物的功能信号 | 定量误差、分层相关；避免大量背景零掩盖活性区误差 |
| 混合估计有效 | 独立高质量目标及低深度/缺失输入 | 实测、序列、融合；同时报覆盖与来源比例 |
| E–G增强联系 | 有方向与功效规则的功能扰动 | 距离、ABC-style、适用gABC、eta=0/1；AP、固定阈值precision/recall、覆盖 |
| 额外组学增益 | 同背景组学和功能标签 | base-only与多组学学习、校准；缺层模式分层 |
| 个体遗传变化 | 同位点跨动物配对、等位测量或变异扰动 | 参考与个体输入、适当基线；位点内变化误差和方向 |
| 育种相关性 | 独立QTL/eQTL、精细定位或功能证据 | 距离/LD/表达/可测性匹配对照；关联支持不等于因果真值 |

不同任务分开设置 `locus_prediction_validated`、`individual_delta_validated`、`link_prediction_validated`。跨位点表现不能证明个体变异方向准确，已有个体基因组研究说明这一风险。[Sasse等](https://pubmed.ncbi.nlm.nih.gov/38036778/)

划分依据研究问题：未见区域用染色体/区块，未见动物用donor/family，未见品种另留品种；不要强行把每个小数据集拆成不可评估的多重五折。训练/调参/校准/测试角色固定，所有预处理只在训练侧拟合。缺类别或有效分组时报告不可评估，不能复制样本补足。

候选内评分比较与端到端发现比较分开。共同候选比较按第7节重算共同分母；端到端还要把预定评价范围内、因候选未发现或不可评分而遗漏的已测试阳性计入漏检。首次atlas-only WGS不宣称新元件发现能力。

precision/recall的运行阈值由训练/校准集冻结；测试集PR曲线可描述，但不能用其反向选择部署阈值。候选发现阈值与功能联系决策阈值分别校准。无适用功能标签就仅报连续分数。eQTL共定位和Hi-C输入不能代替独立功能联系真值。

## 11. 家养动物数据局限与对应处理

| 局限 | 对应处理 | 仍然存在的限制 |
|---|---|---|
| 新动物无表观组 | 同物种/组织序列模型和接触先验 | 不能恢复即时非遗传状态 |
| 整个物种无功能数据 | 明确外推、报告资产缺口 | 无目标验证不宣称适配完成 |
| 组织/阶段/性别覆盖偏倚 | context与输出头限定、分层验证 | 未见状态可能失败 |
| 品种/家系样本少且相关 | 正则化、family/breed划分 | 不能靠更多参数解决样本不足 |
| Hi-C稀疏/参考组织不同 | 单位匹配先验、收缩、近距离策略 | 不能恢复所有个体环结构 |
| TSS注释不完整 | 物理去重、起始证据、明确等权 | 未知TSS仍会遗漏 |
| 参考偏倚、SV/CNV、旁系同源 | callability、映射和影响范围标记 | 未报告变异和复杂剂量仍受限 |
| bulk细胞混合 | 固定估计对象和背景 | DNA不能唯一识别细胞比例 |
| 甲基化低覆盖/平台偏差 | CpG计数、覆盖、E/P分开 | RRBS缺失不等于零甲基化 |
| 技术缺失改变分母 | 三类宇宙ID、共同分母比较 | 条目覆盖不代表真实支持覆盖 |
| 扰动标签稀少/方向混杂 | 标签方向、功效与独立ML | 不能由关联数据冒充因果性能 |

资源局限依具体物种与研究而异，不能说家养动物普遍没有功能数据。已有牛猪鸡图谱有其特定采样范围，也不能称作每个动物九组学完整配套。[Kern等](https://doi.org/10.1038/s41467-021-22100-8)

## 12. 参数与首次发布范围

| 参数/对象 | 初始设置或来源 |
|---|---|
| estimand / target_level | 主输出bulk_proxy；individual/population_mean显式指定 |
| regime | measured/hybrid/genome_only，用户选择并校验实际来源 |
| eta | 0；1是预注册对照 |
| activity_panel | 五种合法组合之一，运行内固定 |
| unit_width / grid_offset | 默认500 bp/0；改变后独立catalog及匹配训练 |
| pi | 可靠起始证据；否则去重可信TSS等权 |
| w、s | w按独立校准拟合；s按训练协议冻结；不猜通用质量权重 |
| r、a_s、gamma_s、d_min | 观测/先验模式及同物种接触拟合，不内置人类最优值 |
| candidate_radius | 启动例5 Mb，须核对物种和验证范围 |
| sequence weights | 实际训练制品，含单位、目标窗口、背景及验证报告 |
| CNN训练起点 | 8192 bp输入、500 bp目标；AdamW lr=1e-3、batch=64、最多50轮、patience=5；不是已调优结果 |
| lambda1/lambda2 | 分组训练侧选择；无标签不拟合 |
| callable/CpG QC | 按实验协议配置并保存；不宣传统一最优阈值 |
| ploidy/phase/callability | 来自样本可靠元数据，不由软件臆造 |
| link threshold | 适用功能标签校准；无标签则缺省 |
| candidate discovery threshold | 仅发现扩展使用，与link threshold分开 |

首次发布实现三模式、标准数据接口、规定常见格式、定量序列适配/参考训练、校准融合、atlas候选、SNV/小indel、组成比较、九组学特征、基础学习/基准和离线示例。全基因组扫描、复杂拷贝重构、任意单倍型组装映射、序列接触训练和复杂后验抽样列研究扩展，未实现时不注册假成功命令。

## 13. 必须成立的手算验收

| 项目 | 输入 | 正确结果 |
|---|---|---|
| 活性 | ATAC=4，H3K27ac=9 | A=6 |
| 多TSS | C=(2,6)，pi=(3/4,1/4) | Cbar=3 |
| 主评分 | A=(4,2,1)，C_G1=(3,2,2)，C_G2=(1,2,6) | G1 eta0=(2/3,2/9,1/9)；eta1=(18/23,4/23,1/23) |
| 融合 | s=1，obs=9，seq=3，w=1/2 | sqrt(40)-1 |
| bulk顺序 | E1两hap=(9,1),(1,9)，E2=(4,4),(4,4)，C全1 | 逐层平均后A=(5,4)，E1 P=5/9 |
| 拷贝研究估计 | 上述相同输入先算每hap支持并相加 | E1 P=3/7；不得替换bulk结果 |
| 共同分母 | 样本a支持(1,1,2)，b支持(1,1,NA) | 原E1分数1/4与1/2不可直接比；共同两单元各1/2；完整Delta=NA |
| 组成变化 | S从(1,1)到(1,2) | E1支持不变，P从1/2到1/3 |
| 同比变化 | S从(1,1)到(2,2) | P不变，总支持翻倍 |
| 甲基化 | 两CpG (m,n)=(1,2),(8,8) | pooled=0.9，site=0.75 |

追加状态验收：全零分母NA；有效零支持为0；技术缺失NA；重复来源不增加评分单元；同一输入改变顺序不改变结果；核心缺失不被ML填补；不支持的SV影响边不继续静默评分。合成验收只证明这些规则的实现正确，不证明预测准确。

<!-- CANONICAL_MODEL_END -->
