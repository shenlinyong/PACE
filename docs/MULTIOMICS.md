# RNA、组蛋白、CTCF 和 DNA 甲基化接口

PACE 区分主公式中的活性/接触与附加组学特征。把不同生物学信号分别保留下来，可以检查它们提供了什么证据，并在有独立功能标签时检验其增益。软件不把所有文件压成一个未经验证的加权和。

## 先登记样本，再决定数据进入哪里

一个项目可以只有 1 个生物学重复，也可以有 2、3 个或更多。每个样本在 `samples.tsv` 中保留独立的 `sample_id`、`donor_id`、`assay`、生物学重复号和技术重复号；`observed_activity.tsv`、`observed_contacts.tsv`、`expression.tsv`、`methylation.tsv` 与 `features.tsv` 使用相同的 ID。不同检测层不要求文件数量相同，但必须说明哪些样本和供体可比较。

主分数只从明确选择的 ATAC/DNase/H3K27ac 活性组合和接触证据构建。RNA-seq、H3K4me1、H3K4me3、H3K27me3、H3K9me3、CTCF、WGBS/RRBS 等可以同时输入，但默认作为命名注释；只有经过独立功能标签和特征契约验证，才进入 ML 分析。不要为了“使用所有数据”把不同量纲直接相乘。

## 支持的数据和角色

| 数据 | 推荐量化对象 | 入口 | 默认处理 |
|---|---|---|---|
| ATAC / DNase | 固定元件窗口 | `observed_activity` | 显式选择的主活性层 |
| H3K27ac | 固定元件窗口 | `observed_activity` | 显式选择的主活性层 |
| H3K4me1 | 元件窗口 | 附加 activity 行或 `features` | 元件状态注释 |
| H3K4me3 | 启动子窗口 | promoter `features`；元件窗口也可留作独立注释 | 启动子状态注释 |
| H3K27me3 / H3K9me3 | 元件/启动子及合适的宽区域 | 附加 activity 行或 `features` | 抑制相关染色质注释 |
| CTCF | 元件、边界、motif 位点 | bigWig / BED adapter / `features` | 结构相关占据与方向注释 |
| WGBS / RRBS | 元件和链特异定义的启动子窗口 | 原始 CpG 计数 `methylation` | 甲基化水平、覆盖及状态 |
| RNA-seq | 基因；有依据时可提供独立启动子使用证据 | `expression` / RNA adapter | 基因 TPM 注释 |
| Hi-C | 元件 anchor 与 TSS 所在的 bin pair | `observed_contacts` | 主接触项 |

H3K27me3 或甲基化并非在所有位置都对应同等抑制效应；CTCF 占据也不能直接换算成接触值。主活性组合只接受 ATAC/DNase/H3K27ac 的支持组合，其余信号作为命名特征进入独立分析。`multiomics.mode: annotate` 是默认模式。

## 1. 追加组蛋白和 CTCF 定量轨道

以下配置与 ATAC/H3K27ac 使用同一个 bigWig adapter：

```yaml
kind: bigwig
track: data/animal1_H3K4me1.bw
units: prepared/catalog/units.tsv
sample_id: animal1_H3K4me1
assay: H3K4me1
unit: normalized_signal
normalization_id: frozen_H3K4me1_protocol
window_id: grid:500:mean
missing_is_measured_zero: false
minimum_callable_fraction: 0.8
```

```bash
PACE prepare --config prepare_h3k4me1.yaml --out prepared/animal1_h3k4me1
```

将输出行合并到 `inputs.observed_activity` 指定的表，同时在 samples 和 sources 中登记真实样本与处理方法。`activity.panel` 不增加 H3K4me1，主分数不会因此换成三层几何平均；该层会出现在 `multiomics_features.tsv.gz`。

对 H3K4me3 等启动子标记，应使用符合研究问题的启动子区间量化，输出为 `features.tsv` 中的 `entity_type=promoter`。不要把 ±2 kb 启动子信号标记成 `grid:500:mean`。需要自定义量化和模型时保留准确的窗口、量纲和处理协议。

## 2. CTCF BED 与 motif 方向

保存为 `prepare_ctcf.yaml`：

```yaml
kind: bed_features
bed: data/ctcf_motifs.bed
units: prepared/catalog/units.tsv
source_id: ctcf_motif_source
evidence_id: ctcf_motif_evidence
feature_prefix: CTCF
entity_type: element
motif_strands: true
```

```bash
PACE prepare --config prepare_ctcf.yaml --out prepared/ctcf
```

`motif_strands: true` 要求 BED 的 strand 列提供 motif 方向。普通 CTCF ChIP 峰通常不携带 motif 方向：只有峰时设置 false，不要拿峰的任意 strand 代替 motif。该接口提取重叠与方向特征，不推断一个经过验证的 loop 概率。

输出 `features.tsv` 指定的 evidence_id/source_id 必须在输入证据和来源表中登记。单条证据可以被多个特征引用，但不能引用不存在的 ID。额外 feature 名称需区分元件、启动子、基因或具体 E–G 边。

## 3. RNA-seq

直接输入的 `expression.tsv` 表头为：

```text
gene_id	sample_id	tpm	status
```

一个示例数据行可写为 `geneA<TAB>animal1_RNA<TAB>12.5<TAB>observed`，其中 geneA 必须存在于启动子/基因目录，sample_id 必须存在于 samples，且 assay 为 RNA。不可用测量保留相应状态和 NA；即使遗留数值没有清空，非 observed 状态也不会成为有效学习特征。

```yaml
inputs:
  expression: prepared/expression.tsv
```

上面是**合入完整 run 配置的片段**，不能单独运行。若已有 transcript TPM，可转换：

```yaml
kind: rna
expression: data/transcript_tpm.tsv
transcript_mapping: prepared/catalog/transcript_mapping.tsv
```

```bash
PACE prepare --config prepare_rna.yaml --out prepared/rna
```

输出的 gene TPM 只按已提供的 transcript-to-gene/TSS 映射汇总，不会根据 gene TPM 反推转录本或 TSS 的使用比例。输出测量与元件预测来自不同个体时，应明确目标是群体注释还是个体分析，不能只靠相同组织名称混合。

## 4. WGBS / RRBS：推荐直接把计数接入主流程

`inputs.methylation` 接收每个 CpG dyad 的计数，不接收百分比 bedGraph 或 summary 表：

```text
chrom	dyad_start0	methylated_count	total_count	sample_id	assay
```

dyad_start0 是正向参考序列 CpG 中 C 的 0-based 位置。每个 dyad/sample/assay 一行；甲基化计数不能超过总计数。相反链同一 CpG 应先正确合并，不能当成两个独立 CpG。Python `merge_stranded_cpg` 提供基于参考序列核验的显式链调用转换。

在完整 run YAML 中加入：

```yaml
inputs:
  methylation: data/cpg_dyad_counts.tsv
methylation:
  minimum_coverage: 5
  promoter_upstream_bp: 2000
  promoter_downstream_bp: 500
  reference_cpg_path: data/reference_cpg_counts.tsv
```

这里 coverage=5 和窗口只是分析协议示例。默认 minimum_coverage 为 1；应根据读深和检测方案预先确定。上游/下游按 TSS 转录方向解释。正链窗口是 `[tss0-upstream, tss0+downstream+1)`，负链交换两侧，起点截断到零。元件窗口来自 units。若末端窗口超出参考，参考 CpG 计数应来自与实际参考边界一致的区间。

`reference_cpg_counts.tsv` 表头：

```text
entity_type	entity_id	n_cpg
```

entity_type 为 `element` 或 `promoter`，entity_id 对应 element_id 或 promoter_id，n_cpg 是**参考序列在同一窗口的总 CpG 数**。旧式 `element_id,n_cpg` 表仍可用于元件。缺少参考 CpG 分母时，coverage_fraction 保持 NA，不把“测到的 CpG 数”误作“所有 CpG 数”。

主流程会自动生成元件和启动子层面的命名特征，并登记证据来源。WGBS 和 RRBS 同时存在时保留 assay 命名空间，不把不同覆盖设计的数据静默平均。RRBS 没有覆盖的区域是未知，不是零甲基化。

| 输出指标 | 含义 |
|---|---|
| `M_site` | 对达到覆盖阈值的 CpG，先求各位点甲基化比例，再等权平均 |
| `M_pooled` | 合并甲基化读数后除以总读数，深度高的 CpG 权重更大 |
| `covered_cpg` | 达到 coverage 阈值的 CpG 数 |
| `cpg_coverage_fraction` | 合格 CpG 数 / 对应参考窗口 CpG 总数 |
| 状态 | 区分真实零甲基化、无 CpG、无覆盖与低覆盖 |

普通亚硫酸氢盐测序无法独立区分 5mC 与 5hmC；解释时需按实验类型说明。

单独查看元件汇总也可以：

```yaml
kind: methylation
counts: data/cpg_dyad_counts.tsv
units: prepared/catalog/units.tsv
minimum_coverage: 5
reference_cpg: data/element_reference_cpg_counts.tsv
```

```bash
PACE prepare --config prepare_methylation.yaml --out prepared/methylation_qc
```

这里的 reference_cpg 是 `element_id,n_cpg` 格式。输出 `methylation_summary.tsv` 用于检查，不要把它再传给要求原始计数的 `inputs.methylation`。需要作为自定义附加特征导入时，显式转换为 features 格式并登记证据。

## 5. 其他表观组学的通用入口

任何已量化、含义明确的附加层可通过 `inputs.features` 导入。必需列：

```text
entity_type	entity_id	feature_name	value	evidence_id	status
```

entity_type 是 element、promoter、gene 或 edge；edge 的 entity_id 使用 `element_id|gene_id`。每个 entity_type/entity_id/feature_name 只能有一条记录。多个重复先按声明规则汇总；不能靠重复行改变特征权重。未知、无效和未测状态保留为缺失，非有限数值不能冒充有效数值。

如需表示 promoter_H3K4me3、methylation_M_site、CTCF_orientation，名称和层级必须与模型训练完全一致。一个通用入口表示软件能保存和使用这个特征，**不表示已经证明它提高家养动物预测准确率**。

## 6. 启用独立 ML 模块

有与目标组织匹配、独立分组的功能扰动标签后，按[训练手册](training.md)训练。真实模型的训练配置需要 `feature_contract`：从适用代表性评分运行的 `ml_feature_contract.json` 读取，并作为 YAML mapping 写入训练配置。它声明活性组合/量纲、接触合同、候选规则、目标层级和附加特征处理。

```yaml
multiomics:
  mode: ml
  model_path: models/liver_classifier
```

这是 run 配置片段。也可独立 `PACE predict-ml --config predict_ml.yaml --out results/ml`。模型训练/推理合同不一致、缺乏可验证范围或无效核心分数时，软件给出明确状态，不输出貌似适用的校准概率。

主 PACE 与附加分类器分开解释，并比较 base-only 和增加表观特征后的独立验证效果。数据使用范围、缺失率和功效不足的负标签同样影响结论。
