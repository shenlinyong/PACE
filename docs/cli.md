# PACE 命令与参数

## 常用入口

```bash
PACE demo --out results/demo
PACE validate --config experiment.yaml
PACE run --config experiment.yaml --out results/experiment
PACE measured --config experiment.yaml --out results/another_run
PACE run --help
```

`run` 和 `measured` 使用同一实测活性模型。`--mode measured` 可保留，但无需选择模式。
所有成功结果写入新的目录；已有目录不覆盖。YAML 路径相对于配置所在目录，CLI 文件路径相对于当前目录。

## 不使用 YAML 的完整示例

下例可从仓库根目录直接运行，数值应与 measured 示例一致：

```bash
PACE run --catalog-dir examples/measured \
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
| --out | 必需的新输出目录 |

## 其他命令

| 命令 | 用途 |
|---|---|
| prepare | BED/GTF、bigWig、cool/mcool、RNA、甲基化和区间特征转标准表 |
| capabilities | 检查背景、接触资产及其验证声明；不会代替完整 validate |
| fit-contact-prior | 用实测接触拟合距离背景，保留真实零和独立测试区域 |
| fit-eta | 用合格功能标签学习可选分配指数 |
| train / predict-ml | 独立的实测特征分类器训练与冻结推断 |
| compare / stability | 在共同背景下比较结果、检查重复一致性 |
| benchmark | 独立功能标签评估和基线比较 |

例如：

```bash
PACE fit-contact-prior --config examples/training/contact.yaml --out models/demo_contact
PACE fit-eta --config examples/measured/config.yaml \
  --eta-labels examples/training/eta_labels.tsv --eta-min-genes 2 --out results/eta_demo
PACE benchmark --config examples/analysis/benchmark.yaml --out results/benchmark_demo
```

这些示例中的数据是合成数据。所有 YAML 默认值见[参数手册](parameters.md)，完整真实数据流程见[中文说明](USER_GUIDE.zh-CN.md)。
