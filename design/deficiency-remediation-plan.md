# Deficiency remediation plan

Status: decisions recorded 2026-11-13. Tracks the deficiencies found while
fixing the test suite, reviewing `R/`, and rewriting the vignettes. Every item
links to a beads issue in `.beads/issues.jsonl`.

## How this relates to the 1.0 plan

`design/fmristore-1.0-plan.md` replaces the frame store with a new engine
(Phases 1-5) and removes the legacy surface (`H5Neuro*`, `H5Parcellated*`,
`LabeledVolumeSet`, latent I/O) in Phase 6. The legacy classes have no real
users outside our own packages, so we do not maintain them (D1). Work splits
into:

- **Green CI and an honest changelog** (Milestone 0).
- **Legacy freeze:** no features, no API cleanup. Fix a legacy bug only if it
  breaks CI or corrupts data for one of our own packages before it migrates.
- **Consumer migration:** move our packages off the legacy surface so Phase 6
  can delete it.
- **1.0 engine requirements:** frame-store durability and integrity findings
  become requirements and storage-law tests of the new engine.

## Decisions

| ID | Decision | Consequence |
|----|----------|-------------|
| D1 | Do not maintain the legacy classes; remove them in 1.0 Phase 6. | Legacy API and docs issues close as won't fix. Legacy correctness bugs are deferred unless a consumer or CI is blocked. |
| D2 | Relax `.lintr`: drop `implicit_integer_linter` and `condition_call_linter`, lint only changed files in CI. | Lint and style gates can pass without a package-wide restyle. |
| D3 | Writers must not overwrite silently. | Applies to the 1.0 engine's writers (default `overwrite = FALSE`). Not retrofitted to legacy writers. |
| D4 | Out-of-mask fill stays `0`. | Moot under D1: only legacy readers fill masks. No `fill =` work. |

## Definition of done (every item)

- [ ] A regression test that fails before the change (for code fixes).
- [ ] Full suite passes against the **installed** package with `NOT_CRAN=true`.
- [ ] `R CMD check` is clean on the hosted matrix.
- [ ] Rd updated (by roxygen at the pinned version).
- [ ] NEWS.md entry for any user-visible change.
- [ ] Beads issue closed with a reason.

## Milestone 0: Trustworthy baseline (do first)

Goal: CI tells the truth, and users know what changed.

- [ ] **Get `main` green again.** [Run 43520115994](https://github.com/bbuchsbaum/fmristore/actions/runs/43520115994) failed on every platform after PR #67 enabled the suite. These were latent failures, not regressions. `fmristore-bfi` (P0)
  - [ ] Ubuntu: fix the malformed fixture at `test-h5_corruption_handling.R:122`. `matrix(rnorm(8 * 5), nrow = 5, ncol = 10)` and the `7 * 5` line under-fill the matrix, which R 4.6 makes an error. The test runs only with `NOT_CRAN=true`.
  - [ ] Windows (11 failures): triage from the check artifact. The visible one is `test-vec-regressions.R:278` (`H5NeuroVec(out)`), likely an open handle blocking a copy or reopen. Related: `fmristore-57h`, `fmristore-71g`.
  - [ ] macOS: triage from the check artifact.
  - [ ] Ubuntu and Windows: the upstream `concat` failure, `fmristore-qsh` (below).
  - [ ] Run local and agent checks with `NOT_CRAN=true`, so that `skip_on_cran()` tests run as they do in CI.
- [ ] **Changelog.** Record the behaviour changes already merged in NEWS.md. `fmristore-2v0` (P1)
- [ ] **Pin dependencies and test the installed package in CI.** Pin `Remotes:` to commits or tags. Add a test convention for comparing neuroim2 objects via `@.Data`. `fmristore-5jx` (P2)
- [ ] **Lint gate.** Add `any::cyclocomp` to `lint.yaml`, drop the two linters per D2, and lint only changed files. Exit: the lint job reports and passes. `fmristore-ky0` (P2)
- [ ] **Style gate.** Check style on changed files only. Exit: the style job passes. `fmristore-ky0`
- [ ] **Upstream `concat` bug.** File it on fmrilatent (or neuroim2's `.concat4D`) with the reproduction, and move the test there. It exercises upstream behaviour only, and fmristore will not own latent I/O after 1.0. Delete it from fmristore in the same change, citing the upstream issue. `fmristore-qsh` (P2)

Exit evidence: a green hosted R CMD check matrix on one pinned dependency set.

## Milestone 1: Legacy freeze

Goal: the legacy surface stays exactly as useful as it is today until it is
deleted, at no maintenance cost.

- [ ] **Announce the freeze.** Mark the legacy classes as superseded in README, the vignettes and their help pages (`lifecycle::badge("superseded")`). Point to the 1.0 plan.
- [ ] **No new legacy work.** Reject features and API changes for `H5Neuro*`, `H5Parcellated*`, `LabeledVolumeSet` and latent I/O.
- [ ] **Deferred correctness bugs.** Fix one only if CI is red because of it, or one of our packages needs that path before it migrates. Otherwise these are deleted with the code:
  `fmristore-cmi` (INT truncation), `fmristore-qt9` (silent overwrite), `fmristore-71g` (constructor modes and leaks), `fmristore-o4x` (label 0), `fmristore-bpi` (`H5NeuroVecSeq` cannot load), `fmristore-e2u` (`LabeledVolumeSet` space), `fmristore-crp` (single timepoint), `fmristore-x7o` (metadata layout), `fmristore-57h` (latent paths, unflushed copy).
- [x] **Closed as won't fix (D1):** argument names `fmristore-nmf`, return values `fmristore-1xb`, `series()` index semantics `fmristore-acy`, chatty output `fmristore-5dm`, summary column names `fmristore-99i`, plain-4D `write_dataset()` `fmristore-xo4`, out-of-mask fill `fmristore-wkp`, `as.matrix`/closed errors `fmristore-rmq`, `series(drop =)` `fmristore-9le`, `x[i]` semantics `fmristore-n2z`, duplicated legacy code `fmristore-skz`.

## Milestone 2: Consumer migration (1.0 Phase 6, started early)

Goal: nothing we own calls the legacy surface, so Phase 6 is a deletion.

- [ ] **Inventory** every call into legacy classes and functions from our packages, with an owner and a replacement for each (1.0 Phase 0 deliverable). Known consumers:
  - `fmridataset`: legacy HDF5 backend;
  - `neurotabs`: parcellated ingestion;
  - `fmrigds`: labeled, latent and parcellated adapters and fixtures;
  - the Python ports.
- [ ] **Read-only decoders.** For each legacy file format we must still read, decide between a read-only decoder (internal or companion package) and a documented refusal.
- [ ] **Migrate** each consumer to the store API, or to the read-only decoders, as the 1.0 engine reaches the needed phase.
- [ ] **Delete** the legacy surface, its tests and vignettes, and the `neuroim2`/`fmrilatent` imports. The deferred M1 issues close with it.

## Milestone 3: 1.0 engine requirements (fold into 1.0 Phases 1-5)

Goal: these become requirements and storage-law tests of the 1.0 engine, not
patches to the RDS-backed store. Map each to its 1.0 section and make sure the
1.0 test suite covers it.

- [ ] **No silent overwrite** (D3): every engine writer defaults to `overwrite = FALSE` and publishes without clobbering.
- [ ] **Writer lock** and lease honoured by rollback, and publication that cannot clobber an existing file (`file.link`). 1.0: Transaction protocol and Concurrency. `fmristore-o0y` (P1)
- [ ] **Durable commit:** fsync of file and directory, and close errors surfaced. 1.0: Atomicity and Process death. `fmristore-b62` (P2)
- [ ] **Fingerprint checks.** Validate shape, dtype and fingerprint on `source_open()`, using sub-second mtime or a content token. 1.0 Phase 1 (hardened `h5_array_source`). `fmristore-0t2` (P2)
- [ ] **Manifest encoding.** Replace `serialize()`/`unserialize()` manifests with a canonical, schema'd encoding. 1.0: Identity and integrity. `fmristore-2un` (P2)
- [ ] **Dtype contract:** float16/bfloat16 sizing, NA in float32, and native read types. 1.0: Dtype and missing-value contract. `fmristore-lr4` (P3)
- [ ] **Recovery provenance.** Prove that a crashed, never-committed store directory is ours (marker file). 1.0: Recovery and garbage collection. `fmristore-css` (P3)
- [ ] **O(1) shard appends** (catalog index instead of opening every shard). 1.0: catalog. `fmristore-0h2` (P3)

Interim rule: if any of these bites a user before 1.0, fix it in the current
store with a regression test, and keep the 1.0 requirement open.

## Milestone 4: Readability and documentation hygiene (new code only)

Apply to the 1.0 engine and the code that survives Phase 6. Do not spend it on
legacy files.

- [ ] **Runnable examples** for every surviving export (tempfile-based, no `\dontrun{}`). `fmristore-6z4` (P3)
- [ ] **One condition style:** the 1.0 structured condition hierarchy (classed `cli_abort()`), enforced by lint. `fmristore-fub` (P3)
- [ ] **Roxygen.** Pin the roxygen2 version in CI and regenerate `man/` once. Several Rd files were edited by hand in this pass. `fmristore-5jx`

## Suggested order and size

| Order | Work | Size |
|-------|------|------|
| 0 | Green `main` (`fmristore-bfi`): fixture fix, Windows and macOS triage | S-M |
| 1 | Rest of M0: NEWS, CI pinning, lint/style gates, move the `concat` test upstream | S-M |
| 2 | M1 freeze notices | S |
| 3 | M2 consumer inventory | M |
| 4 | M3 inside the 1.0 phases; M2 migrations as the engine allows | (1.0 schedule) |
| 5 | M2 deletion (1.0 Phase 6) | M |

M4 runs alongside steps 4-5.
