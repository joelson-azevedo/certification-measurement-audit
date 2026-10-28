# Reproduction tiers and exact local commands

**Synthetic portfolio candidate with approved author and license decisions.** The included project-specific synthetic delivery is licensed CC BY 4.0; external publication of this exact candidate still awaits separate authorization. Do not use the research outputs for operational certification decisions. [E01] [E08]

## Environment

The observed interpreter is CPython 3.12.14 on Windows x64. Numerical dependencies are NumPy 2.3.5, pandas 3.0.1 and SciPy 1.18.1; tests use pytest 9.1.1. PDF and image inspection use pypdf 6.10.0 and Pillow 12.3.0. Exact constraints are in [quick requirements](environment/requirements-quick.txt); dependency versions are engineering environment observations, not scientific result estimates.

No dependency download, network request or service is needed by the verification program itself. A Python interpreter and dependencies must already exist or be installed first. The candidate does not contain Python, wheels, compiled third-party libraries or a wheelhouse. Consequently, a new offline machine cannot install the requirements from this package alone. Obtain approved dependencies separately; an offline installation can use an existing wheelhouse. No tools are assumed to be preinstalled silently.

### PowerShell

From this candidate's root, create and populate an environment. The installation below is deliberately offline and requires a separately prepared wheelhouse:

```powershell
python --version
python -m venv .venv
.venv/Scripts/python.exe -m pip install --no-index --find-links $wheelhouse -r environment/requirements-quick.txt
.venv/Scripts/python.exe -I -X utf8 analysis/verify_release.py --run-tests
```

If a compatible environment is already installed, use its interpreter for the final command. `$wheelhouse` is a user-supplied directory of approved wheels; no such directory is embedded or downloaded here. The release verifier refuses missing/mismatched scientific dependencies instead of substituting another backend.

### POSIX shell

The code uses relative paths and portable Python APIs. These equivalent commands are provided; POSIX execution has not been tested in this Windows validation:

```sh
python3 --version
python3 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links "$WHEELHOUSE" -r environment/requirements-quick.txt
.venv/bin/python -I analysis/verify_release.py --run-tests
```

`WHEELHOUSE` must refer to approved wheels appropriate to that interpreter and platform. Exact dependency pins do not guarantee byte-identical floating-point behavior across platforms. No fixed runtime promise is made; actual local timings are reported by the verifier.

## Tier 1: quick verification

Run the command above from this candidate, not the original audit directory. It performs these checks:

1. Verify the curated package manifest and the original delivered-data manifest before analysis; check path containment and complete file coverage.
2. Validate the delivered schemas, membership, response options and chronological histories; independently rebuild scores and decisions from the proposed data.
3. Compare counts, raw/decision discrepancies and observed first-attempt rates to retained original results. Compare all corrected item-rest correlations, alpha and KR-20 point values using the preserved functions. Alpha bootstrap intervals are retained evidence, not recalculated by this quick path. [MP02]
4. Recompute the primary observed-score correspondence, fixed-reference conditional agreement and integrated original-rule agreement using frozen item parameters and the saved provisional link. This verifies calculations conditional on those estimates; it does not refit the latent model, establish invariance or recreate the uncertainty sample. [MP08] [MP09]
5. Resolve every indexed claim to its packaged source cell/configuration/excerpt, verify its formatting and context, validate links, figures, PDF copies and source-module provenance.
6. Import the selected numerical modules and run the packaged tests, including independent controlled fixtures and release checks. Fixture success does not establish validity of the examination model.

Results are written only to the untracked `verification/` directory; pytest uses contained `.test_runs/` fixtures. These runtime files are excluded from the static manifest. Original proposed inputs, retained expected results and published documents are never overwritten. Package hash failure stops verification. The retained manifest establishes internal consistency, not external authentication.

The verification is intentionally narrower than full estimation. It checks retained DIF and bootstrap evidence by source identity/hash and cell semantics; it does not refit every regional model, repeat parallel analysis, run a complete residual simulation or recalibrate the paired bootstrap. The [claim registry](docs/claim_registry.json) and [results guide](results/README.md) distinguish retained evidence from regenerated quick outputs.

## Tier 2: complete analytical replay is unavailable as distributed

This package contains unchanged core numerical implementations and exact recorded configurations, but not the original stage executors, full processed releases, historical stage inventories or complete parameter/link uncertainty caches. Those dependencies are integral to the original stage-gated pipeline; pretending that their absence is a minor path change would overstate reproducibility.

The [full replay plan](docs/FULL_REPLAY_PLAN.md) documents stage order, required inputs/outputs, computationally intensive procedures and seeds. It is a controlled-workspace reference plan, **not a runnable full-pipeline command in this distribution**. Releasing the missing components or engineering and independently validating a portable replacement workflow would require additional work and another review. Tier 1 success does not close that gap.

The author explicitly approved publication with the Tier 2 limitation. It is a disclosed scope limit, not an unresolved human-approval decision.

The approved executive PDFs are copied, not rebuilt. Regenerating their exact bytes would additionally require the original production source, embedded-font-compatible licensed system fonts, ReportLab and Poppler. Those publication-production dependencies are not required for reading or verifying the retained PDFs and are not silently installed.

## Clean local validation boundary

The candidate is tested in a separate copied directory using a fresh isolated venv. Necessary installed dependency distributions are copied locally into that venv; no network is used. Package code/data imports and paths must resolve inside the copied candidate; interpreter and dependency paths are recorded internally. This tests independence from the original audit files. It does not establish fresh wheel installation, POSIX behavior or reproducibility on unspecified hardware. See [publication status](docs/RIGHTS_STATUS.md) for the separate exact-candidate authorization still required before external distribution.

[E01]: docs/claims.md#e01
[E08]: docs/claims.md#e08
[MP02]: docs/claims.md#mp02
[MP08]: docs/claims.md#mp08
[MP09]: docs/claims.md#mp09

## Metadata and network controls

The verifier also parses the complete citation, checks approved author/affiliation and mixed-license scope, compares original MIT legal wording, checks the exact approved AI disclosure and confirms the license inventory covers all static files. A formal CFF schema validator is unavailable offline; parsing and explicit schema-oriented checks are documented without claiming formal certification. Python socket/DNS operations are denied by an audit hook. The selected pytest suite runs in the same process so that network control remains active. This is a Python-level safeguard, not an operating-system sandbox for arbitrary native code. No verifier component is permitted to start a subprocess.
