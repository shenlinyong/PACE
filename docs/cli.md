# PACE 命令与参数

## 常用入口

```bash
pace demo --out results/demo
pace validate --config experiment.yaml
pace run --config experiment.yaml --out results/experiment
pace measured --config experiment.yaml --out results/another_run
pace run --help
```

`run` 和 `measured` 使用同一实测活性模型。`--mode measured` 可保留，但无需选择模式。
默认写入新目录。`--force` 先把已有 PACE 结果保存到相邻备份目录，再发布新结果；不允许替换任意输入目录。YAML 路径相对于配置所在目录，CLI 文件路径相对于当前目录。

## 不使用 YAML 的完整示例

下例可从仓库根目录直接运行，数值应与 measured 示例一致：

```bash
pace run --catalog-dir examples/measured \
  --species synthetic --assembly toy_assembly --tissue toy_tissue \
  --profile demonstration --target-level individual \
  --panel ATAC H3K27ac --contact-mode observed --contact-scale toy_contact \
  --no-include-promoters --out results/direct_measured
```

真实数据请更换目录、背景和测量尺度，使用 `--profile research`。
`--catalog-dir` 自动读取存在的标准表；同一文件的显式参数优先。该目录应包含 units、promoters、candidates、samples、sources 及活动/接触表。

| 参数 | 含义 |
|---|---|
| --config | 完整运行 YAML，命令行参数可覆盖对应项 |
| --activity / --observed-activity | 实测信号表，可同时包含多个检测层及重复 |
| --contacts / --observed-contacts | 实测元件—TSS 接触表 |
| --resolved-activity / --resolved-contacts | 有完整来源记录的已解析结果；活性仅接受实测或实测聚合 |
| --species / --assembly / --tissue | 必须与样本及先验范围一致的生物学背景 |
| --target-level | individual 单供体；population_mean 多供体等权汇总 |
| --panel | ATAC、DNase、H3K27ac 单层，或 ATAC H3K27ac、DNase H3K27ac |
| --contact-mode | observed 默认；prior_only 显式距离先验；shrinkage 显式接触收缩 |
| --contact-prior | 适用的已拟合接触先验；使用先验或收缩时必需 |
| --contact-resolution / --contact-scale | 统一分辨率和信号尺度 |
| --contact-normalization / --contact-balancing / --contact-window | 完整接触测量定义 |
| --contact-reliability / --reliability-source | 接触收缩的实测权重 [0,1] 及其依据 |
| --allow-prior-fallback | 明确允许缺失接触回退到已提供的先验 |
| --eta | auto 默认；无合格功能标签使用零，也可预先固定 [0,1] 内数值做敏感性分析 |
| --eta-labels / --eta-model | 自动学习所需标签，或已冻结且适用的校准资产 |
| --expression / --methylation / --features | RNA、CpG 计数和其他表观注释入口 |
| --ml-model | 已独立训练的附加分类器，输出与主分数分列 |
| --out | 输出目录 |
| --force | 保留旧结果备份后替换已识别的 PACE 输出 |
| --pseudocount | auto、none 或 powerlaw；只能用同尺度先验正则化实测接触 |
| --partial-policy | withhold 默认；conditional 显式兼容条件主分数 |
| --support-bounds | 有来源说明的缺失最终支持上下界 |
| --prior-preset | abc_human；仅 research + prior_only，目标背景未验证 |

## 其他命令

| 命令 | 用途 |
|---|---|
| init | 生成简短配置和输入空表 |
| prepare-pairs | 从候选目录直接生成 E–TSS 查询表 |
| merge-tables | 合并样本规范表并检查重复 |
| normalize-activity | 原始窗口计数按文库量和窗口长度做 CPM 密度归一化 |
| prepare-promoter-weights | 从实测启动子信号生成固定 pi |
| fit-prior | 直接从 cool/mcool 拟合先验，或通过 --config 使用原表格接口 |
| prepare | BED/GTF、bigWig、cool/mcool、RNA、甲基化和区间特征转标准表 |
| capabilities | 检查背景、接触资产及其验证声明；不会代替完整 validate |
| fit-contact-prior | 用实测接触拟合距离背景，保留真实零和独立测试区域 |
| fit-eta | 用合格功能标签学习可选分配指数 |
| train / predict-ml | 独立的实测特征分类器训练与冻结推断 |
| compare / stability | 在共同背景下比较结果、检查重复一致性 |
| benchmark | 独立功能标签评估和基线比较 |

例如：

```bash
pace fit-contact-prior --config examples/training/contact.yaml --out models/demo_contact
pace fit-eta --config examples/measured/config.yaml \
  --eta-labels examples/training/eta_labels.tsv --eta-min-genes 2 --out results/eta_demo
pace benchmark --config examples/analysis/benchmark.yaml --out results/benchmark_demo
```

这些示例中的数据是合成数据。所有 YAML 默认值见[参数手册](parameters.md)，完整真实数据流程见[中文说明](USER_GUIDE.zh-CN.md)。

## Optional sparse-data controls

`pace-livestock` is a lowercase alias of `pace`. No uppercase executable is installed.

| Option | Meaning |
|---|---|
| --allow-cross-context-prior | Explicitly permit a same-species, same-assembly prior from another tissue |
| --minimum-tss-weight | Filter TSSs by original weight, using one set per gene |
| --missing-tss-policy strict or drop_missing | Keep all positive-weight TSSs, or remove those missing any planned contact |
| --minimum-retained-tss-weight | Require this much original weight after filtering; default 0.9 |

Set assay-specific offsets in YAML with `activity.pseudocounts`; they default to zero.
For `fit-prior`, omitted `--min-distance` uses the matrix resolution.
