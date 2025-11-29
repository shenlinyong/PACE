# PACE: Mathematical Framework and Methods

## 1. The ABC Model Foundation

The Activity-by-Contact (ABC) model predicts enhancer-gene regulation based on:

```
                    A(E) × C(E,G)
ABC Score(E,G) = ─────────────────────
                 Σ[A(e) × C(e,G)]
                 e∈5Mb
```

Where:
- **E**: Candidate enhancer element
- **G**: Target gene
- **A(E)**: Activity of enhancer E
- **C(E,G)**: Contact frequency between E and G promoter
- **5Mb**: All candidate elements within 5 Mb of gene G

## 2. PACE Enhanced Activity Calculation

### 2.1 Multi-layer Activity Score

PACE extends activity calculation to integrate multiple data types:

```
A(E) = Aggregation(S₁, S₂, ..., Sₙ) × (1 - I(E))
```

Where:
- **Sᵢ**: Individual signal (ATAC, H3K27ac, H3K4me1, etc.)
- **Aggregation**: Configurable function
- **I(E)**: Inhibitory score from repressive signals

### 2.2 Aggregation Methods

#### Geometric Mean (Original ABC style):
```
A = (∏ᵢ Sᵢ)^(1/n)
```

#### Weighted Geometric Mean (PACE default):
```
A = ∏ᵢ (Sᵢ^wᵢ)
```

#### Weighted Sum:
```
A = Σᵢ (wᵢ × Sᵢ)
```

### 2.3 Inhibitory Score

Combines repressive marks and DNA methylation:

```
I(E) = α × Methylation(E) + β × H3K27me3(E) + γ × H3K9me3(E)
```

### 2.4 Complete Activity Formula

```
A(E) = [∏ⱼ (Activating_Signal_j)^wⱼ] × [1 - Σₖ (wₖ × Inhibitory_Signal_k)]
```

**Example:**
```
A(E) = (ATAC^1.5 × H3K27ac^1.0 × H3K4me1^0.8) × (1 - 0.5×Methylation)
```

## 3. Contact Estimation

### 3.1 Hi-C Based (if available):
```
C(E,G) = Hi-C_KR_normalized(bin_E, bin_G)
```

### 3.2 Power-law Distance Decay (default):
```
C(E,G) = Scale / (Distance + Pseudocount)^γ

Parameters:
- Scale = 5.959
- γ = 1.024
- Pseudocount = 5000 bp
```

## 4. Expression Integration

### 4.1 Expression Filter

Only consider genes with expression ≥ threshold:
```
Gene_set = {G : Expression(G) ≥ θ}
Default θ = 1 TPM
```

### 4.2 Expression Weight Methods

```
Binary:   W(G) = 1 if Expression(G) ≥ θ else 0
Linear:   W(G) = Expression(G) / max(Expression)
Log:      W(G) = log₂(Expression(G) + 1) / max(log values)
```

## 5. Final PACE Score

```
                    A(E) × C(E,G) × W_expr(G)
PACE Score(E,G) = ─────────────────────────────
                   Σ[A(e) × C(e,G) × W_expr(G)]
```

## 6. Signal Weight Recommendations

| Signal | Type | Default Weight | Biological Basis |
|--------|------|----------------|------------------|
| ATAC/DNase | Accessibility | 1.5 | Open chromatin, TF accessible |
| H3K27ac | Active enhancer | 1.0 | Active enhancers/promoters |
| H3K4me1 | Enhancer mark | 0.8 | Poised + active enhancers |
| H3K4me3 | Promoter mark | 0.5 | Active promoters |
| H3K36me3 | Transcription | 0.3 | Transcribed gene bodies |
| H3K9ac | Active chromatin | 0.6 | Generally active regions |
| H3K27me3 | Repressive | 0.5 | Polycomb repression |
| H3K9me3 | Heterochromatin | 0.5 | Constitutive heterochromatin |
| Methylation | Inhibitory | 0.5 | Transcriptional silencing |
| TF binding | Activating | 0.3 | Regulatory proteins |

## 7. Workflow Overview

```
Input Data                    Processing                     Output
───────────                   ──────────                     ──────

ATAC-seq ─────┐
              │
H3K27ac ──────┼──► Peak Calling ──► Candidate Regions
              │         │
H3K4me1 ──────┤         │
              │         ▼
H3K4me3 ──────┼──► Activity Calculation ──┐
              │                           │
Methylation ──┤                           │
              │                           ▼
TF ChIP ──────┘                    E-G Pair Generation
                                          │
Hi-C ─────────────────────────────────────┤
                                          │
                                          ▼
RNA-seq ──────────────────────► Expression Filter/Weight
                                          │
                                          ▼
eQTL ─────────────────────────► Validation (optional)
                                          │
                                          ▼
                                   PACE Predictions
                                          │
                                          ▼
                              ┌───────────────────────┐
                              │ • Enhancer coordinates │
                              │ • Target genes         │
                              │ • PACE scores          │
                              │ • Distance metrics     │
                              │ • Confidence levels    │
                              └───────────────────────┘
```

## 8. Usage Scenarios

| Scenario | Configuration | Data Required |
|----------|--------------|---------------|
| **Minimal** | ATAC only | ATAC-seq |
| **Standard** | ATAC + H3K27ac | + H3K27ac ChIP |
| **Enhanced** | Multi-histone | + H3K4me1/me3 |
| **Expression-aware** | + RNA-seq | + RNA-seq |
| **Repression-aware** | + Methylation | + WGBS/RRBS |
| **Mechanism-focused** | + TF binding | + TF ChIP-seq |
| **Genetically-validated** | + eQTL | + eQTL data |
| **Comprehensive** | All available | All data types |

## 9. Parameter Reference

### Distance Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Max E-G distance | 5,000,000 bp | Maximum enhancer-gene distance |
| Power-law γ | 1.024 | Distance decay exponent |
| Power-law scale | 5.959 | Distance decay scale factor |
| Pseudocount | 5,000 bp | Avoid division by zero |

### Filtering Thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| ABC/PACE score | 0.02 | Prediction confidence threshold |
| Min expression | 1 TPM | Gene expression filter |
| eQTL p-value | 1e-5 | eQTL significance |
| Peak p-value | 0.1 | MACS2 peak calling |

### Candidate Region Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Peak extend | 250 bp | Extension from peak summit |
| Region size | 500 bp | Total candidate region size |
| N strongest peaks | 150,000 | Maximum peaks to retain |
