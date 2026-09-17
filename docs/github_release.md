# Maintainer release checklist

PACE source is maintained at [shenlinyong/PACE](https://github.com/shenlinyong/PACE).
A source commit, package version, GitHub release, container publication and software
DOI are separate records. Refer to an artifact as published only after it actually
exists at the stated location.

Before submitting a release commit:

1. Run the regression suite with the required optional dependencies, the public
   documentation checker, lint and distribution checks.
2. Run all three example modes from the installed wheel; check Conda and Docker
   jobs for that commit.
3. Review the diff for unintended files, generated results, private data and large
   model assets. Include only redistributable source, documents and fixtures.
4. Record the tested environment and known limitations in the release notes.
5. Push without rewriting existing history, then inspect hosted CI and the published
   documentation. Tag a release only after its acceptance checks pass.

```bash
git status --short
git diff --check
python scripts/check_public_docs.py
```

Pin source commit and dependencies for a scientific analysis. Real biological
model assets need their own model cards, validation data scope and redistribution
licenses; software regression success cannot replace those checks.
