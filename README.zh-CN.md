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

[完整公式说明](docs/```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```.md)给出**全部展开的通用总公式、三种模式分别展开的总公式、每个符号和数值算例**。缺失是 NA，真实零是 0；分数是相对支持，不能解释为因果概率或基因表达贡献比例。

## 其他表观组学与结果

软件支持 RNA-seq、ATAC/DNase、H3K27ac、H3K4me1、H3K4me3、H3K27me3、H3K9me3、CTCF、WGBS/RRBS 和 Hi-C。[多组学手册](docs/MULTIOMICS.md)列出接口、角色、甲基化覆盖和启动子窗口、CTCF motif 方向及 ML 接入。附加层默认保留为可追溯注释，不随意乘入主分数。

首先查看 `gene_summary.tsv` 的候选覆盖、`qc_report.json` 的缺失原因，再看 `scores.tsv.gz` 的分数。`partial` 代表只对可评分子集归一化。`resolved_activity.tsv`、`resolved_contacts.tsv` 和证据目录说明每个值来自哪里。`run_manifest.json` 与 `resolved_config.yaml` 用于重现结果。

多物种参考版本、表观数据不齐、Hi-C 分辨率差异、多 TSS、供体重复及鸡的性染色体倍性等都有明确处理规则。这些设计使框架便于适配家养动物；具体物种、组织和品种的预测效果仍需独立实验验证。本仓库**没有附带已验证的真实家养动物序列权重**。

[完整中文说明书](docs/USER_GUIDE.zh-CN.md) · [全部参数](docs/parameters.md) ·
[输入字段](docs/data_dictionary.md) · [训练](docs/training.md) ·
[差异比较](docs/comparison.md) · [局限](docs/limitations.md)。
