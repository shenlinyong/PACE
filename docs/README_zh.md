# PACE（Prediction of Activity-based regulatory Connections for Enhancers）

**利用增强子活性、染色质接触和启动子注释，预测增强子与候选靶基因之间的调控联系，并明确报告证据缺失。**

作者与维护者：**shenlinyong — 申林用（Linyong Shen），西北农林科技大学**。

[English](../README.md) · [安装说明](INSTALLATION.md) · [逐步教程](TUTORIAL.md) · [参数与依据](PARAMETERS.md) · [结果字段](IO_FORMATS.md)

## 软件解决什么问题

对一个物种、一个基因组版本和一种组织或细胞类型，软件根据增强子活性与启动子接触计算候选增强子–基因联系的相对支持。输入可以是已定量的数值表，也可以是候选区间、注释和基因组信号文件。主模型无需训练标签或 GPU。

家养动物研究中，不同组织的测定类型、测序覆盖和启动子注释往往不一致，组织匹配的 Hi-C 也可能缺少。因此，软件提供缺失值处理、多 TSS 汇总、接触来源标记和独立证据状态，帮助研究者在信息有限时保留候选，并识别需要进一步验证的环节。使用其他物种时必须提供相应参考文件；没有“所有动物通用的最优参数”。

输出用于排序和实验候选筛选。**分数不是调控概率，也不是已校准的假发现率。**

## 直接与原始 ABC 比较

比较对象为 Fulco 等人在 2019 年发表的 [Activity-by-Contact 模型](https://doi.org/10.1038/s41588-019-0538-0)及其 [NG2019 代码](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019)。原始 ABC 已提供平均接触和距离模型替代方案，不能把“没有匹配 Hi-C 仍可运行”单独列为新贡献。

| 改动 | 为什么改 | 对家养动物数据的用途 | 使用边界 |
| --- | --- | --- | --- |
| 缺失感知的平移几何活性 | 区分未测量与测得零，保留可用测定 | 可声明单模态方案或处理局部覆盖缺失 | 不能推算不存在的测量；低信号假阳性需评估 |
| 显式测定尺度与默认等权 | 避免不同信号单位和未经验证的权重混在一起 | 使不同组织、数据预算下的分析方案可记录 | 等权不是从家养动物数据学得的最优值 |
| 接触可靠性加权 | 将合格的观测/期望接触与距离先验结合 | 可以保留稀疏接触或借用组织的信息与来源 | 期望曲线、质量和组织匹配关系需独立提供 |
| 去重后的多 TSS 加权平均 | 避免仅因转录本多而放大接触 | 保留已注释的可变启动子 | 不会恢复漏注释的启动子或基因 |
| 增强子侧候选靶基因分配 | 同时考虑一个增强子附近的多个候选基因 | 把靶基因竞争背景纳入评分 | 基因目录缺失也会影响分配结果 |
| 残余支持量接口 | 明确表达候选集合外支持的独立估计 | 允许记录有根据的遗漏支持 | 默认未知，计算时取零；软件不自动估计 |
| 分数与证据质量分开 | 高相对分数不代表所有输入都充分 | 可同时检查优先级与证据缺口 | 证据状态不是功能真阳性标签 |

靶基因分配和多 TSS 思想引用 [Hecker 等（2023），generalized ABC / STARE](https://doi.org/10.1093/bioinformatics/btad062)，其中本软件采用 TSS 加权平均实现。具体公式与限制见[模型比较](ABC_COMPARISON.md)和[公式说明](FORMULA.md)。

## 参数为何这样设置

| 参数 | 默认值 | 依据与调整方式 |
| --- | --- | --- |
| 活性模式 | `missing_geometric` | 区分零与缺失；采用 `log1p` 后的加权平均 |
| ATAC/DNase 与 H3K27ac 权重 | 各 1 | 作为等权设计先验；在活性 JSON 中用 `weight` 明确指定 |
| 信号尺度 | 示例为 1 | 生物学数据需按预先规定的参考规则选择 `scale`，不能直接把不同测定原始计数当作同单位 |
| 靶分配指数 | 1 | 使用靶分配；`--competition-power 0` 用于移除这一组件的消融 |
| 接触可靠性 | 缺少条件时用于混合的值为 0 | 回退到距离先验，同时保留原因；已知质量用 `contact_reliability` 提供 |
| TSS 权重 | 去重后均分 | 无独立启动子使用数据时的假设；有数据可提供 `tss_weight` |
| 顺式窗口 | 严格小于 5 Mb | 宽范围候选搜索起点，不是动物特异最优距离 |
| 距离指数、尺度与偏移 | 约 1.02424、5.95945、5,000 bp | 可追溯的计算起点，未宣称为独立拟合的家养动物参数 |
| 抑制强度 | 0 | 主分析不强制使用甲基化或抑制性标记；显式扩展通过核心 API |
| 证据阈值 | 0.5 | 操作性输入质量规则，不是统计显著性 |
| 输出分数阈值 | 0.02 | 示例筛选值；新物种或组织需要独立校准 |
| 峰扩展 | 峰顶两侧各 250 bp | 形成约 500 bp 的初始候选；候选策略必须一致 |
| MACS2 有效基因组大小 | 通用配置中为占位值 | 必须按实际物种、组装与比对方案修改，鸡不能直接套用哺乳动物值 |

[完整参数表](PARAMETERS.md)逐项列出实际入口。表格命令行、单样本 YAML 和 Snakemake 可调范围不同；有些保留的 YAML 字段未被转交给当前预测入口，不能只增加一个字段就认为计算已改变。

## 使用 Conda 安装

先安装 Linux 环境、Bash、Git 与 Conda，例如 [Miniforge](https://github.com/conda-forge/miniforge#install)。Python 和软件依赖由环境文件一次性安装。

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate pace
```

此环境包含 Python 3.11、NumPy、pandas、PyYAML、SciPy、matplotlib、pyBigWig、bedtools、samtools 和 pytest。无需逐个预装这些包。

已有候选峰时不需要 MACS2；直接运行逐步命令时不需要 Snakemake。若需要从已比对 reads 开始自动调用 MACS2，请另建工作流环境：

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f workflow/envs/pace-env.yml
conda activate pace-workflow
```

`.hic` 和 `.cool` 分别需要可选的 `hic-straw` 和 `cooler`；BEDPE 不需要额外读取器。参见[依赖用途与安装故障处理](INSTALLATION.md)。软件不执行 FASTQ 比对或原始 Hi-C 矩阵构建。

## 第一个可运行例子

在项目根目录、`pace` 环境中执行：

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/quickstart/predictions.tsv
```

预期生成 **6 条增强子–基因联系**；其中 2 条因活性缺失而保留为空值。输入是 9 条增强子–TSS 记录，多 TSS 汇总后行数减少是正常结果。

如需直接从随附的合成 reads 运行完整示例：

```bash
bash example/run_example_direct.sh
```

此命令完成候选构建、信号定量、评分、筛选和 QC，输出位于 `example/results/Example_Sample/`，完整表含 **12,000 条候选联系**。合成示例验证软件行为，不证明生物学性能。

验证安装：

```bash
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

## 处理自己的动物组织数据

准备同一组装版本的 FASTA/GTF/染色体长度、可及性峰和经过预处理的 ATAC/DNase 信号；有同组织 H3K27ac 时可加入。用 `prepare_tss.py` 生成去重的 TSS 表，用 `pace_neighborhoods.py` 定量，然后用 `pace_predict.py` 预测。每一步的完整命令、YAML/样本表模板及可选 Hi-C/RNA 参数都在[逐步教程](TUTORIAL.md)。

要求所有文件染色体名称和坐标一致；`NA` 表示缺失，`0` 表示测得零。跨物种或跨组装坐标不会自动转换。原始 reads/kb 不等于跨文库深度归一化；需要固定尺度和局部 QC 时，优先准备数值表使用 `scripts/pace.py`。

## 如何解释结果

首先保存全部候选结果，再另存阈值筛选表。重点查看：

| 字段 | 含义 |
| --- | --- |
| `PACE.Score` | 当前候选背景内的相对支持 |
| `contact_gene` / `target_share` / `raw_support` | 分数的组成部分 |
| `TargetGeneTSSs` | 参与汇总的 TSS 坐标 |
| `contact_state` | 距离先验、匹配组织、借用组织或缺少质量信息等状态 |
| `evidence_status` / `evidence_reasons` | 输入充分性及具体限制 |
| `unscored_candidates` | 同基因所给候选中无法评分的数量 |

仅有距离先验或缺少质量信息时，出现 `provisional` 通常是预期行为；不要为了得到“充分”标签而将未知质量填为 1。`sufficient_input_evidence` 也不等于实验验证阳性。所有列和伴随文件见[结果说明](IO_FORMATS.md)。

引用时记录仓库及实际 Git 提交，同时按方法来源引用 ABC 和 generalized ABC/STARE。当前仓库未指定 PACE 论文 DOI。问题请提交到 [GitHub Issues](https://github.com/shenlinyong/PACE/issues)，附命令、环境和可复现的小输入。
