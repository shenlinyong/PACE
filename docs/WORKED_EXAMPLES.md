# Independent worked calculations

These small analytical examples test the current formula; they are not biological
performance results. The numerical tests cite these same hand-derived quantities.

## Fixed-panel activity and TSS contact

ATAC=4 and H3K27ac=9 give activity sqrt(4*9)=6. If H3K27ac is measured zero,
activity is zero. If it is unavailable, activity is NA under the two-assay panel.
Contacts (2,6) with promoter weights (3/4,1/4) give Cbar=3. A missing required TSS
contact does not permit reweighting the remaining TSS to weight one.

## Endpoint and continuous allocation

Let activities be (4,2,1), G1 contacts (3,2,2), and G2 contacts (1,2,6).
Then B for G1 is (3/4,1/2,1/4), and B for G2 is (1/4,1/2,3/4).

| Exponent | G1 support | G1 PACE | G2 PACE |
|---|---|---|---|
| 0 | (12,4,2) | (2/3,2/9,1/9) | (2/7,2/7,3/7) |
| 1 | (9,2,1/2) | (18/23,4/23,1/23) | (2/15,4/15,9/15) |

For eta=1/2, G1 support is (6sqrt(3),2sqrt(2),1) and G2 support is
(2,2sqrt(2),3sqrt(3)). Divide each vector by its own sum. Every listed denominator
is exactly the sum of those supports. A completely zero vector produces NA scores.

## Bulk aggregation and composition

Two equally weighted measured replicates of E1 have assay pairs (9,1) and (1,9). Two equally weighted measured replicates of E2 each have (4,4).
The bulk calculation first averages assays, producing activities 5 and 4. With
identical contact and eta zero, E1 receives 5/9. Computing support separately per
replicate and summing gives 3/7 instead; that is a different, unimplemented estimand.

Supports (1,1,2) and (1,1,NA) give original E1 shares 1/4 and 1/2. On the common
first two elements, both recomputed shares are 1/2; a complete-background difference
remains unavailable. Changing supports (1,1) to (1,2) reduces the first share from
1/2 to 1/3 even though its own support is unchanged.

## Continuous calibration check

Two equally weighted pairwise margins eta-1/2 and 1/2-eta yield the convex loss
`[softplus(1/2-eta)+softplus(eta-1/2)]/2`. Its unique minimum is eta=1/2 and its value
is log(2). This independently checks an interior solution, not just the two endpoints.

[Scoring tests](../tests/test_canonical_core.py) ·
[Continuous calibration tests](../tests/test_canonical_allocation.py) ·
[Formula](FORMULA.md) · [Run a tutorial](TUTORIAL.md).
