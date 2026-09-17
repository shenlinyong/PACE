# Repository delivery and future releases

The authorized target is the existing public repository `https://github.com/shenlinyong/PACE`.
Development preserves Git history while keeping one current scientific implementation. The installable canonical package is
versioned independently from the public name PACE. Uploading code is not a formal GitHub release,
PyPI publication or DOI registration.

Before pushing, run the documented tests and wheel/demo checks, review the staged diff, and
verify that only code, documentation and small synthetic fixtures are included. The repository
ignores generated results, private input directories, local environments and large model formats.

```bash
git status --short
git diff --check
git push origin main
```

No force push or history rewrite is needed. Existing branches are retained. A future release
should pin its source commit and runtime environment, rerun CI, review real model cards and
licenses separately, then create an authorized tag/archive. The package name `pace-livestock`
was not present at the PyPI JSON endpoint during this development session (2026-09-17); that
check does not reserve the name or publish the package.

The source tree, source distribution and wheel use the same `pace_livestock` implementation,
including its self-contained demo generator. Retired implementations are available through
Git history only. Check public documentation with `python scripts/check_public_docs.py` before
publication, and verify the published default-branch files after pushing.
