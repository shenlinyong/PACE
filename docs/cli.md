# PACE 命令行使用手册

安装后直接使用 `PACE --mode measured|hybrid|genome 参数`。也可以写成
`PACE measured 参数`、`PACE hybrid 参数`、`PACE genome 参数`。
`PACE run --help` 列出全部文件参数；`PACE --help` 列出准备、训练和分析命令。
命令别名 `pace` 和原来的 `pace-livestock` 使用同一套实现。

## 安装

需要 Python 3.11 或更新版本。在克隆的仓库中执行：

```bash
./install.sh --prefix "$HOME/.local" --python python3
export PATH="$HOME/.local/bin:$PATH"
PACE --version
```

脚本将隔离环境安装到 `PREFIX/share/pace/venv`，将命令放入 `PREFIX/bin/PACE`。
系统 Python 过旧时，使用 `--python /path/to/python3.11` 或更新解释器。
它拒绝覆盖非本脚本管理的环境和已有命令。安装需要可用的包索引或本地依赖缓存。

读取 bigWig、cool/mcool、BCF 加 `--extras io`；真实 CNN 训练和推理加
`--extras sequence`，也可以使用 `--extras io,sequence,ml`。PyTorch 的 CPU/GPU
发行包应与计算环境匹配。基础安装可以运行全部合成示例，不下载模型权重。

已有 Python/Conda 环境的用户可直接运行 `python -m pip install .`，或安装
`python -m pip install '.[io,sequence,ml]'`。安装后在该环境中运行 `PACE`。
无需在每次使用时输入 Python 脚本路径。

## 三种模式的真实数据命令

下面的 `sheep`、`my_assembly`、`liver` 和文件路径是需要替换的项目示例，
必须与输入表和模型清单逐字匹配。`prepared/catalog` 包含 `units.tsv`、
`promoters.tsv`、`candidates.tsv`，以及适用的 `sources.tsv`、`evidence.tsv`。
这些文件定义非重叠评分单元、去重启动子、固定候选关系及证据来源。
格式见[数据字典](data_dictionary.md)，从 BED/GTF、bigWig、cool/mcool
生成它们的方法见[输入准备](input_preparation.md)。

实测组学模式：

```bash
PACE --mode measured \
  --catalog-dir prepared/catalog \
  --samples data/samples.tsv \
  --activity data/observed_activity.tsv \
  --contacts data/observed_contacts.tsv \
  --species sheep --assembly my_assembly --tissue liver \
  --panel ATAC H3K27ac --contact-scale depth_normalized_contact \
  --out results/measured
```

活动和接触参数接收规范 TSV/TSV.GZ 表，不直接接收 FASTQ、BAM 或 bigWig。
真实轨道先通过 `PACE prepare --config preparation.yaml --out prepared/track`
提取到匹配窗口。`--sources`、`--evidence` 可以覆盖目录中的同名表。
单层 ATAC、DNase 或 H3K27ac 使用相应 `--panel`，不会静默填补另一层。

混合模式：

```bash
PACE --mode hybrid \
  --catalog-dir prepared/catalog \
  --samples data/samples.tsv \
  --activity data/observed_activity.tsv \
  --contacts data/observed_contacts.tsv \
  --reference data/reference.fa --sequence-model models/sequence \
  --fusion-model models/fusion \
  --species sheep --assembly my_assembly --tissue liver \
  --out results/hybrid
```

`--fusion-model` 是独立目标数据校准过的融合资产。没有适用融合资产时，按
既有来源选择规则解析证据，不任意设置实测/预测权重。模型单位、窗口和组织
必须匹配。需要个体化序列时，加下面基因组命令中的 VCF、可调用区域与倍性参数；
实测供体须与该个体一致。

基因组预测模式：

```bash
PACE --mode genome \
  --catalog-dir prepared/catalog \
  --reference data/reference.fa \
  --vcf data/animal.vcf.gz --sample-id animal1 --individual-id animal1 \
  --callable data/callable.bed --ploidy data/ploidy.tsv \
  --sequence-model models/sequence --contact-prior models/contact \
  --species sheep --assembly my_assembly --tissue liver \
  --out results/genome
```

`genome` 是 `genome_only` 的命令行别名；无配置文件时自动选择 `prior_only`
接触模式。catalog 目录中不能混入实测活动/接触表。仅做参考基因组背景预测时，
去掉 VCF、sample/individual、callable 和 ploidy 参数。FASTA、候选注释与适用的
定量模型及接触先验仍然必需。仅有 FASTA 不能产生经过验证的组织特异预测。

## 自动估计和固定 eta

默认是 `--eta auto`：无适用功能标签时实际使用 **0**。提供具备方向、阴性功效、
背景与数据划分信息的规范标签表后自动估计连续参数：

```bash
PACE measured --config project.yaml \
  --eta-labels functional_labels.tsv --out results/calibrated
```

`PACE fit-eta` 接受相同参数，额外要求提供功能标签。输出包含评分和可重复使用的
`eta_calibration.json`。以后在适用范围一致的分析中冻结它：

```bash
PACE measured --config another_animal.yaml \
  --eta-model results/calibrated/eta_calibration.json --out results/frozen
```

也可以预先指定 `--eta 0.35`。手动数值与校准输入互斥。关于估计目标、标签格式、
最少信息量、排除原因和测试集隔离，见[eta 校准说明](eta_calibration.md)。

## 配置、路径与结果

`--config project.yaml` 可省略；直接参数与配置并用时，直接参数优先。
命令行路径相对于当前工作目录，YAML 内路径相对于 YAML 所在目录。
`--catalog-dir` 自动读取存在的标准同名 TSV/TSV.GZ 表，显式文件参数再覆盖；
同一表同时存在压缩和未压缩版本会报歧义错误。

主结果是 `scores.tsv.gz`，另有逐基因分母、证据表、QC、运行清单、实际配置和
eta 校准记录。实际使用的 eta 写入终端摘要和运行清单；`auto` 不是一个数值结果。
退出码 0 表示命令成功，输入错误等返回 2。已有输出目录不会覆盖。成功但部分
候选不可评分时，QC 和每条记录会说明原因；技术缺失保持 NA。

## 可以立即运行的合成检查

在仓库根目录完成安装后：

```bash
PACE measured --config examples/measured/config.yaml --out results/check_measured
PACE hybrid --config examples/hybrid/config.yaml --out results/check_hybrid
PACE genome --config examples/genome_only/config.yaml --out results/check_genome
PACE fit-eta --config examples/measured/config.yaml \
  --eta-labels examples/training/eta_labels.tsv --eta-min-genes 2 \
  --out results/check_eta
```

最后一个例子只有两个合成基因，因此显式使用 2；真实分析默认要求至少三个有信息基因；
还需足够独立分组、每折训练信息和稳定的留出增益。这个小示例应返回零回退，
不能用来证明非零参数获得了验证。仓库不附带已验证的真实家养动物权重。

## 显式测量与注释参数

`PACE run --help` 是当前参数的直接查询入口。新增 CLI 与 YAML 对应关系：

| CLI | YAML |
|---|---|
| `--chrom-sizes` | catalog.chrom_sizes_path |
| `--candidate-radius` | catalog.candidate_radius_bp |
| `--contact-resolution` | contact.resolution |
| `--contact-normalization` | contact.normalization_id |
| `--contact-balancing` | contact.balancing |
| `--contact-window` | contact.window_id |
| `--reference-cpg` | methylation.reference_cpg_path |
| `--methylation-min-coverage` | methylation.minimum_coverage |
| `--promoter-upstream` / `--promoter-downstream` | methylation.promoter_upstream_bp / promoter_downstream_bp |
| `--eta-min-groups` | allocation.minimum_groups |
| `--eta-validation-folds` | allocation.validation_folds |

`--catalog-dir` 会读取 prepare 导出的 run_catalog_config.yaml 中的目录参数；
显式 CLI 参数随后覆盖。完整真实配置和逐项解释见[中文手册](USER_GUIDE.zh-CN.md)。
