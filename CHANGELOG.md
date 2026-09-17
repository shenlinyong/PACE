# PACE changes

## 17 September 2026

- Use PACE consistently as the software and model name; the mathematical score is italic PACE(E,G).
- Replace the separate prediction implementations with the validated shared scoring kernel.
- Preserve observed zero, missing measurements and unknown input quality as separate states.
- Add distinct-TSS contact aggregation, enhancer-target allocation and explicit contact provenance.
- Correct the earlier expression-weighted formula in all retained documentation: normalize activity × contact first, then multiply by the gene expression weight. The weight does not enter the denominator.
- Keep RNA as context in current primary predictions; archive incompatible earlier ML and analysis tools.
- Supply a quantified-table example, command-line smoke workflow, regression tests and current input documentation.
- Update the bundled read example to the supported missing-aware activity mode.

The release preserves the scoring equations used by the current manuscript analysis. The earlier expression-weighted model is documented for comparison and is not silently reintroduced into primary predictions.
