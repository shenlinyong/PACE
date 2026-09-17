# PACE changes

## 17 September 2026

- Clarify PACE as a general framework particularly suited to incomplete livestock datasets; expand all six original-ABC comparisons with parameter rationale and worked calculations.
- Replace the overview diagram with formula comparisons and executable input scenarios; expand Conda prerequisites and the Chinese manual.

- Rebuild the user manual around the original ABC comparison, parameter rationale, Conda installation, executable examples and output interpretation; add a Chinese guide.
- Provide separate direct-workflow and Snakemake Conda recipes, exact Linux package locks and a Conda-only setup helper.
- Match filter/QC log wildcards to their outputs so the documented Snakemake workflow can build and execute its DAG.

- Use PACE consistently as the software and model name; the mathematical score is italic PACE(E,G).
- Replace the separate prediction implementations with the validated shared scoring kernel.
- Preserve observed zero, missing measurements and unknown input quality as separate states.
- Add distinct-TSS contact aggregation, enhancer-target allocation and explicit contact provenance.
- Correct the earlier expression-weighted formula in all retained documentation: normalize activity × contact first, then multiply by the gene expression weight. The weight does not enter the denominator.
- Apply the same correction to the archived sensitivity-analysis calculator, with regression checks for non-unit expression weights, optional expression input and zero scores. Expression-weighted legacy recalculations can therefore differ from earlier outputs.
- Keep RNA as context in current primary predictions; archive incompatible earlier ML and analysis tools.
- Supply a quantified-table example, command-line smoke workflow, regression tests and current input documentation.
- Update the bundled read example to the supported missing-aware activity mode.

The release preserves the scoring equations used by the current manuscript analysis. The earlier expression-weighted model is documented for comparison and is not silently reintroduced into primary predictions.
