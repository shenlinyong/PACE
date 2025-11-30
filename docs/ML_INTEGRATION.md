# PACE Machine Learning Integration

## Overview

PACE includes an experimental machine learning module that can enhance enhancer-gene predictions by learning optimal feature combinations from validated data.

## Model Selection

### Default Model: Gradient Boosting

PACE默认使用 **Gradient Boosting Classifier** (梯度提升分类器)，原因如下：

| 特性 | 说明 |
|------|------|
| **可解释性** | 提供特征重要性排序，帮助理解哪些表观遗传信号最重要 |
| **非线性关系** | 可以捕捉特征之间的复杂非线性交互（如ATAC×H3K27ac的协同效应） |
| **处理不平衡数据** | E-G预测中正样本远少于负样本，GB对此有较好的鲁棒性 |
| **不需要特征缩放** | 决策树基础的模型对特征尺度不敏感 |
| **防止过拟合** | 通过学习率、树深度等参数控制模型复杂度 |
| **文献验证** | 在类似基因组学问题中被广泛使用并验证有效 |

### 备选模型: Random Forest

也可以选择 **Random Forest Classifier**：

```bash
python scripts/pace_ml.py train \
    --model_type random_forest \
    ...
```

| 特性 | Gradient Boosting | Random Forest |
|------|-------------------|---------------|
| 训练速度 | 较慢（顺序构建） | 较快（并行构建） |
| 预测精度 | 通常更高 | 略低 |
| 过拟合风险 | 需要调参 | 更鲁棒 |
| 特征重要性 | ✓ | ✓ |

### 为什么不用深度学习？

对于E-G预测问题，我们选择传统ML模型而非深度学习：

1. **样本量限制**: 验证的E-G对通常只有几百到几千个，不足以训练深度网络
2. **可解释性需求**: 研究者需要理解哪些特征驱动预测
3. **计算资源**: 传统ML模型不需要GPU，更易部署
4. **性能足够**: 在表格型数据上，GB/RF的性能通常不逊于神经网络

### 相关文献中的模型选择

| 工具 | 模型 | 参考 |
|------|------|------|
| TargetFinder | Random Forest | Whalen et al. 2016 |
| JEME | Elastic Net | Cao et al. 2017 |
| EP-DNN | Deep Neural Network | Li et al. 2019 |
| ABC Model | 无ML（公式驱动） | Fulco et al. 2019 |
| **PACE** | Gradient Boosting | 本工具 |

## Motivation

The standard PACE/ABC model uses a fixed formula:
```
Score = (Activity × Contact) / Σ(Activity × Contact)
```

While this works well in general, the optimal combination of features may vary by:
- Tissue type
- Species
- Available data types
- Biological context

The ML module learns these optimal combinations from validated E-G pairs.

## How It Works

### 1. Feature Engineering

The ML module extracts features from PACE predictions:

| Feature Category | Examples |
|-----------------|----------|
| **Epigenomic signals** | ATAC, H3K27ac, H3K4me1, H3K4me3 |
| **Distance features** | log_distance, inverse_distance, is_proximal |
| **Contact features** | contact, activity_contact_product |
| **Expression features** | log_expression, is_expressed |
| **Derived features** | H3K27ac/H3K4me1 ratio, activity_sum |

### 2. Model Training

Uses Gradient Boosting or Random Forest classifiers:
- Cross-validated to prevent overfitting
- Class balancing for imbalanced data
- Automatic feature importance ranking

### 3. Prediction

Outputs:
- `ML.Score`: ML-predicted probability (0-1)
- `Combined.Score`: Weighted average of ABC.Score and ML.Score

## Usage

### Prerequisites

```bash
pip install scikit-learn
```

### Training a Model

```bash
python scripts/pace_ml.py train \
    --predictions results/Sample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --validation data/eqtl_validated_pairs.tsv \
    --output models/my_model.pkl \
    --model_type gradient_boosting \
    --n_estimators 100 \
    --max_depth 5 \
    --balance_classes
```

#### Validation Data Format

TSV file with validated E-G pairs:

```tsv
enhancer	gene	validated
chr1:1000-1500	BRCA1	1
chr1:5000-5500	TP53	1
chr2:3000-3500	MYC	0
```

- `enhancer`: Region in `chr:start-end` or `chr_start_end` format
- `gene`: Gene symbol matching TargetGene column
- `validated`: 1 for validated positive, 0 for negative (optional)

#### Validation Data Sources

| Source | Description |
|--------|-------------|
| **eQTL** | SNPs affecting gene expression |
| **CRISPRi** | CRISPR interference perturbation |
| **Hi-C loops** | Chromatin loop anchors |
| **Capture-C** | Targeted chromatin capture |
| **Published benchmarks** | ENCODE, Roadmap validated sets |

### Applying a Model

```bash
python scripts/pace_ml.py predict \
    --predictions results/NewSample/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
    --model models/my_model.pkl \
    --output results/NewSample/Predictions/EnhancerPredictions_ML.tsv \
    --abc_weight 0.5
```

#### `--abc_weight` Parameter

Controls how to combine ABC and ML scores:
- `0.0`: Use only ML score
- `0.5`: Equal weight (default)
- `1.0`: Use only ABC score

## Model Configuration

### In config.yaml

```yaml
ml_integration:
  enabled: false                    # Enable ML in workflow
  model_type: "gradient_boosting"   # Model type
  use_pretrained: false             # Use pre-trained model
```

### Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--model_type` | gradient_boosting | Model algorithm |
| `--n_estimators` | 100 | Number of trees |
| `--max_depth` | 5 | Maximum tree depth |
| `--learning_rate` | 0.1 | Learning rate (GB only) |
| `--cv_folds` | 5 | Cross-validation folds |
| `--balance_classes` | False | Downsample negatives |

## Output Interpretation

### Training Output

```
Training Results
==================================================
CV AUC: 0.847 ± 0.023
Samples: 50000
Features: 15

Top Feature Importance:
  contact: 0.2341
  log_distance: 0.1892
  H3K27ac: 0.1456
  ATAC: 0.1234
  activity_contact_product: 0.0987
```

### Feature Importance

Indicates which features are most predictive:
- Higher importance = more predictive
- Use to understand tissue-specific regulation

## Best Practices

1. **Sufficient validation data**: Need at least 50-100 validated pairs
2. **Class balance**: Use `--balance_classes` if positives << negatives
3. **Cross-validation**: Always check CV AUC for model quality
4. **Feature selection**: Start with default, then refine based on importance
5. **Species transfer**: Models trained on one species may transfer to related species

## Limitations

- Requires validated E-G pairs for training
- May overfit on small datasets
- Not integrated into main Snakemake workflow (standalone tool)
- Experimental feature - use with caution

## Transfer Learning

For livestock species with limited validation data:

1. Train on human/mouse data
2. Fine-tune on available livestock data

```python
from ml_integration import PACEMLPredictor, TransferLearning

# Train on human
human_model = train_pace_ml_model(human_data, ...)

# Transfer to pig
transfer = TransferLearning(human_model)
transfer.fine_tune(pig_X, pig_y)

# Predict
probs, preds = transfer.predict(new_pig_X)
```

## Citation

If using the ML module, please cite both PACE and scikit-learn.
