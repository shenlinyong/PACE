# PACE

PACE 使用**实测调控活性**和来源明确的启动子接触证据，计算家养动物候选调控元件对基因的相对支持。

[完整中文使用说明](docs/USER_GUIDE.zh-CN.md) · [公式逐项解读](docs/FORMULA.zh-CN.md) ·
[数据字典](docs/data_dictionary.md) · [参数](docs/parameters.md) · [English](README.md)

## 下载与配置环境

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
pace --version
```

没有 Conda 时，可用 Python >=3.11：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install '.[io,ml]'
```

Docker 安装和运行：

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work pace:local \
  demo --out /work/results/docker_demo
```

无需 GPU。wget 下载、Docker 路径和独立目录安装见[安装说明](docs/INSTALLATION.md)。

## 先运行一个完整示例

```bash
pace demo --out results/demo
pace validate --config examples/measured/config.yaml
pace run --config examples/measured/config.yaml --out results/measured
pace measured --help
```

仓库示例为合成数据，可直接运行，用于检查安装和熟悉流程。真实分析需要替换为自己的实测数据，设置 `execution_profile: research`。
默认输出目录须尚不存在；`--force` 保留旧 PACE 结果备份后替换。`pace run` 和 `pace measured` 是同一个实测活性入口。

## 需要准备什么

| 数据 | 要求 |
|---|---|
| 活性 | 至少有 ATAC、DNase 或 H3K27ac 中一个合格实测层；固定单层或受支持的双层组合 |
| 接触 | 同背景且可比较的 Hi-C 等接触数据；缺少 Hi-C 时可显式提供适用距离先验并注明来源 |
| 候选目录 | 固定评分单元、去重 TSS、候选元件—基因关系 |
| 样本信息 | 真实供体、实验类型、生物/技术重复、物种、组装和组织 |
| 额外组学 | RNA、H3K4me1/3、抑制性标记、CTCF、WGBS/RRBS 为可选注释或独立分类器特征 |

缺失活性保留 NA。只有 H3K27ac 时，可明确选择 `--panel H3K27ac`；已选择双层后不能按元件临时删层。
实际文件格式和从 BED/GTF、bigWig、cool/mcool 到标准表的命令见[输入准备](docs/input_preparation.md)。

## 总公式

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

A_star 为实测活性的几何平均；Cbar 为多启动子综合接触；B 为可选的跨基因分配项。
没有合格功能验证数据时使用 eta=0。主分数要求完整计划候选支持；缺失时另列条件分数和敏感性区间。
[完整展开公式](docs/FORMULA.zh-CN.md)说明每一个符号、重复处理和接触策略。

## 分析真实实验

```bash
pace run --config experiment.yaml --out results/experiment
```

[完整中文手册](docs/USER_GUIDE.zh-CN.md)提供可修改的 YAML、逐项参数、真实文件组织和结果解释。
YAML 路径相对于配置文件；CLI 路径相对于当前目录。单只动物使用 `individual`，多供体汇总需要显式设置 `population_mean`。

重点查看 `scores.tsv.gz`、`gene_summary.tsv`、`qc_report.json` 和 `run_manifest.json`。
PACE 是相对支持分数，不能当成因果概率；软件测试通过也不能代替家养动物中的独立功能验证。

[多组学接口](docs/MULTIOMICS.md) · [验证与比较](docs/comparison.md) · [能力边界](docs/limitations.md) ·
[论文范围调整](docs/MANUSCRIPT_SCOPE.zh-CN.md) · [迁移说明](docs/migration.md)

维护者：申林用（Linyong Shen）。[MIT 许可](LICENSE)。

## 稀疏实验数据与便捷准备

支持同尺度幂律伪计数、近对角校正、缺失分母敏感性区间、原始区域汇总，以及启动子实测信号权重。`pace init` 生成简短配置；`prepare-pairs` 和 `merge-tables` 代替手写连接/合并脚本；`fit-prior` 直接拟合 cool/mcool。

[实际操作](docs/PRACTICAL_WORKFLOW.md) · [意见核查与调整](docs/METHOD_REVIEW.zh-CN.md) · [独立验证方案](docs/ROBUSTNESS_VALIDATION.md)。
