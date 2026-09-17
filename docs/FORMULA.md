# Scoring equations

PACE (Prediction of Activity-based regulatory Connections for Enhancers) computes relative enhancer–gene support. [Parameters](PARAMETERS.md) gives the exact defaults and actual configuration interfaces; [Notation](NOTATION.md) maps symbols to software fields.

## 1. Enhancer activity

For assay $i$, use a nonnegative preprocessed signal $S(E,i)$ and fixed positive scale $a_i$, giving $x(E,i)=S(E,i)/a_i$. The observation indicator $m(E,i)$ is one for a finite measurement and zero for missing input. With assay prior $w_i$ and local quality $q(E,i)$,

$$
A(E)=\left\{\exp\left[\frac{\sum_i m(E,i)q(E,i)w_i\log(1+x(E,i))}{\sum_i m(E,i)q(E,i)w_i}\right]-1\right\}\exp[-\kappa I(E)].
$$

Measured zero contributes to the mean; missing input does not. With no effective observations, activity is undefined. Unknown quality uses weight 1 for provisional activity but stays unknown in quality output. Scales must be chosen independently of evaluation labels and held fixed across planned sensitivity runs. The primary setting is $\kappa=0$; the core API's optional inhibitory extension needs explicitly scaled [0,1] inputs.

This differs from the unshifted geometric activity used in the original ABC profile. With scaled signals 0 and 8 and equal quality/weights, the shifted activity is 2, whereas the unshifted geometric mean is 0. This numerical change is a design choice requiring evaluation, not evidence of biological superiority.

## 2. Contact with explicit reliability

Every enhancer–TSS pair requires a finite positive prior $P(d)$. File adapters use

$$
P(d)=\mathrm{Scale}(d+\mathrm{Pseudocount})^{-\gamma}.
$$

With a measured contact $H(E,t)$, compatible positive expected contact $D(d(E,t))$, and independently assessed reliability $\lambda(E,t)$,

$$
C_{\mathrm{adj}}(E,t)=P(d(E,t))\left[1-\lambda(E,t)+\lambda(E,t)\frac{H(E,t)}{D(d(E,t))}\right].
$$

The observation is usable only when these quantities are valid and its source is declared `matched` or `surrogate`. Otherwise the mixing reliability is zero and the limitation is reported. At reliability zero the prior is recovered; at one, qualified observed/expected contact fully modulates the prior. A measured zero can reduce contact. A missing edge cannot be replaced by zero.

Expected contact and reliability are inputs, not quantities automatically fitted by PACE. Their sample, resolution, normalization and source tissue must be recorded.

## 3. Integrate distinct promoters

For distinct TSSs $\mathcal T(G)$ with prespecified use weights,

$$
C(E,G)=\sum_{t\in\mathcal T(G)}\pi(G,t)C_{\mathrm{adj}}(E,t),
\qquad \sum_{t\in\mathcal T(G)}\pi(G,t)=1.
$$

Uniform use is the default when independent promoter information is unavailable. Weights are set over the whole available gene catalogue before distance filtering; remaining weights are not renormalized separately for each enhancer. Duplicate transcripts sharing a TSS do not add promoter weight. Missing annotation cannot be recovered by this formula.

## 4. Allocate enhancer support across targets

For the candidate genes $\mathcal G(E)$ of one enhancer,

$$
B(E,G)=\frac{C(E,G)}{\sum_{g\in\mathcal G(E)}C(E,g)},
\qquad R(E,G)=A(E)C(E,G)B(E,G)^\eta.
$$

The default $\eta=1$ uses allocation; $\eta=0$ removes it for a component ablation. Missing candidate genes can affect $B$. Allocation and the multiple-TSS principle are attributed to [Hecker et al. (2023)](https://doi.org/10.1093/bioinformatics/btad062); the weighted-average implementation above is the one used here.

When every contact for an enhancer is zero, target allocation is undefined. A finite activity with known zero contact still gives zero raw support; missing activity remains undefined.

## 5. Normalize within genes

$$
\mathit{PACE}(E,G)=\frac{R(E,G)}{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}R(e,G)+U(G)}.
$$

$\mathcal E^{\mathrm{obs}}(G)$ contains supplied candidates with finite raw support. Other supplied candidates remain in the output and contribute to `unscored_candidates`. No finite support, or a zero total denominator, gives an undefined score.

$U(G)$ is optional independent residual support in the same units as $R$. There is no implemented estimator for it. Missing $U$ is set to zero only for computation and flagged `unassigned_mass_unknown`; an explicit zero is a declared assumption. The default therefore scores the observed candidate background rather than correcting incomplete enhancer discovery.

With one TSS, $\eta=0$, $U=0$ and exactly the same supplied $A$ and $C$, this is the original ABC **normalization rule**. It does not make the full preprocessing pipelines equivalent. RNA expression is not a multiplier in this score.

## 6. Report evidence separately

$$
Q(E,G)=\min\{Q_A(E),Q_C(E,G),Q_T(G),Q_{\mathrm{cat}}(G)\}.
$$

For declared assay layers, activity quality is the assay-prior-weighted quality of observed measurements divided by total planned prior weight. Missing planned measurements contribute zero quality; unknown quality on a present layer makes the result unknown. Contact quality is the TSS-weighted mixing reliability. TSS quality uses the conservative minimum across alternatives and retains unknowns. Catalogue quality must be independently supplied.

At the default evidence threshold 0.5, all four components must be known and meet the threshold, with no unscored supplied gene candidates, for `sufficient_input_evidence`. An undefined score or zero activity/TSS quality gives `insufficient`; other cases are `provisional`. A quality designation does not establish a functional true positive.

The core's optional resampling summary reports empirical score quantiles and the fraction of independent reruns with a score. Missing edges count in that fraction's denominator. These are stability summaries, not posterior probabilities or causal confidence intervals.
