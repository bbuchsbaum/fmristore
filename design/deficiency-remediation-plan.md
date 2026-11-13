# Deficiency remediation plan

Status: proposed, 2026-11-13. Tracks the deficiencies found while fixing the
test suite, reviewing `R/`, and rewriting the vignettes. Every item links to a
beads issue in `.beads/issues.jsonl`.

## How this relates to the 1.0 plan

`design/fmristore-1.0-plan.md` replaces the frame store with a new engine
(Phases 1-5) and removes the legacy surface (`H5Neuro*`, `H5Parcellated*`,
`LabeledVolumeSet`, latent I/O) in Phase 6. This plan therefore sorts work
three ways:

- **Fix now (0.x):** correctness and data-safety bugs in code people use
  today. Small, tested, shipped in patch releases.
- **Fold into 1.0:** frame-store durability and integrity work that the 1.0
  engine must solve anyway. Do not patch the RDS-backed store beyond
  regressions.
- **Decide before investing:** API redesign of the legacy surface. Worth it
  only if that surface outlives Phase 6 (decision D1).

## Decisions needed

| ID | Question | Recommended default |
|----|----------|---------------------|
| D1 | Will the legacy classes survive 1.0 (or live on in a companion package)? | Freeze their API. Fix bugs and document, no renames. Revisit if a companion package is chosen. |
| D2 | Lint policy: relax `.lintr` or restyle the package? | Relax: drop `implicit_integer_linter` and `condition_call_linter`, lint only changed files in CI, restyle opportunistically. |
| D3 | Change the default of writers that overwrite silently? | Add `overwrite`, default `TRUE` with a deprecation warning for one release, then `FALSE`. |
| D4 | Out-of-mask voxels: keep `0` or offer `NA`? | Keep `0` as the default, add `fill =`, document both. |

## Definition of done (every item)

- [ ] A regression test that fails before the change.
- [ ] Full suite passes against the **installed** package (not only `load_all()`).
- [ ] `R CMD check` is clean apart from the known upstream failure (`fmristore-qsh`).
- [ ] Rd updated (by roxygen at the pinned version, or by hand to match).
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
- [ ] **Lint gate.** Add `any::cyclocomp` to `lint.yaml` and apply D2. Exit: the lint job reports and passes. `fmristore-ky0` (P2)
- [ ] **Style gate.** Run `styler` on changed files only (or as decided in D2). Exit: the style job passes. `fmristore-ky0`
- [ ] **Upstream `concat` bug.** Report it to fmrilatent/neuroim2 and link the upstream issue. Remove the failing test only once upstream is fixed or the behaviour is explicitly declared unsupported. `fmristore-qsh` (P2)

Exit evidence: a green hosted R CMD check matrix on one pinned dependency set.

## Milestone 1: Legacy correctness and data safety (0.x patch)

Goal: no silent wrong answers or lost data in the classes people use now.

Silent data loss or corruption:

- [ ] **Integer truncation.** `as_h5(data_type = "INT")` must reject values that don't fit, and `?as_h5` must document the float32 default. `fmristore-cmi` (P2)
- [ ] **Silent overwrite.** Add `overwrite` to `as_h5()`, `to_nih5_vec()`, `write_vec()` and `write_labeled_vec()`, per D3. `fmristore-qt9` (P2)
- [ ] **Constructors.** Open readers with `mode = "r"`, validate before opening, and close on error. `write_labeled_vec()` validates before truncating. `fmristore-71g` (P2)
- [ ] **Cluster label 0.** Either reject 0 at write time or treat it as background everywhere, with one ID set. `fmristore-o4x` (P2)

Broken or wrong results:

- [ ] **`H5NeuroVecSeq`** loads files written by `neurovecseq_to_h5()`, and `detect_h5_type()` maps them. `fmristore-bpi` (P1)
- [ ] **`LabeledVolumeSet` space.** Set the space slot so `dim()` and `series()` work. `fmristore-e2u` (P2)
- [ ] **Single-timepoint runs** work in `write_cluster_result()` and `summarize_by_clusters()`, and `summary_fun` is not forced to take `na.rm`. `fmristore-crp` (P2)
- [ ] **Metadata layout.** `write_parcellated_scan_h5()` writes cluster metadata where `H5ParcellatedMultiScan()` reads it, and the summary writer checks `n_time`. `fmristore-x7o` (P2)
- [ ] **Latent file names and copies.** `write_dataset(LatentNeuroVec, "x.h5")` honours the path, and `write_dataset()` flushes before copying. `fmristore-57h` (P2)
- [ ] **`x[i]` on `H5ParcellatedScan`.** One meaning for a single subscript (mask indices), with an error when out of range. `fmristore-n2z` (P2)
- [ ] **Summary column names.** Store `cluster_ids` on every summary dataset and name columns `cluster_<id>` everywhere. `fmristore-99i` (P3)
- [ ] **`series(drop = FALSE)`** is honoured on `H5NeuroVec` (the `H5NeuroVol[, , ]` half is already fixed). `fmristore-9le` (P3)
- [ ] **`as.matrix(H5NeuroVec)`** works, and reads after `close()` give a clear "object is closed" error. `fmristore-rmq` (P3)

Documentation of semantics (cheap, high value):

- [ ] **Out-of-mask fill.** Document that out-of-mask voxels read as `0`; add `fill =` per D4. `fmristore-wkp` (P3)
- [ ] **`series(x, i)` indices.** Document mask-index vs volume-index semantics in `?series`. `fmristore-acy` (P3)

## Milestone 2: Frame store hardening (fold into 1.0 Phases 1-5)

Goal: these become requirements and storage-law tests of the 1.0 engine, not
patches to the RDS-backed store. Map each to its 1.0 section and make sure the
1.0 test suite covers it.

- [ ] **Writer lock** and lease honoured by rollback, and publication that cannot clobber an existing file (`file.link`). 1.0: Transaction protocol and Concurrency. `fmristore-o0y` (P1)
- [ ] **Durable commit:** fsync of file and directory, and close errors surfaced. 1.0: Atomicity and Process death. `fmristore-b62` (P2)
- [ ] **Fingerprint checks.** Validate shape, dtype and fingerprint on `source_open()`, using sub-second mtime or a content token. 1.0 Phase 1 (hardened `h5_array_source`). `fmristore-0t2` (P2)
- [ ] **Manifest encoding.** Replace `serialize()`/`unserialize()` manifests with a canonical, schema'd encoding. 1.0: Identity and integrity. `fmristore-2un` (P2)
- [ ] **Dtype contract:** float16/bfloat16 sizing, NA in float32, and native read types. 1.0: Dtype and missing-value contract. `fmristore-lr4` (P3)
- [ ] **Recovery provenance.** Prove that a crashed, never-committed store directory is ours (marker file). 1.0: Recovery and garbage collection. `fmristore-css` (P3)
- [ ] **O(1) shard appends** (catalog index instead of opening every shard). 1.0: catalog. `fmristore-0h2` (P3)

Interim rule: if any of these bites a user before 1.0, fix it in the current
store with a regression test, and keep the 1.0 requirement open.

## Milestone 3: API consistency (gated on D1)

If D1 = freeze (recommended), do only the boxes marked *freeze*. Otherwise do
all of them, each with a `lifecycle` deprecation cycle.

- [ ] *freeze:* **Argument names.** Document the argument names per function in one reference table. Full option: canonical `file` / `compression` / `chunk_dims` with deprecated aliases. `fmristore-nmf` (P2)
- [ ] *freeze:* **Writer return values.** Document what each writer returns and that a returned H5File must be closed with `$close()`. Full option: writers return the path invisibly. `fmristore-1xb` (P2)
- [ ] **Plain 4D via `write_dataset()`.** `write_dataset(<NeuroVec>)` without `clusters` writes an `H5NeuroVec` file, or the docs name `as_h5()` as the writer. `fmristore-xo4` (P3)
- [ ] **Quiet defaults.** Normal reads and writes print nothing unless `verbose = TRUE`. `fmristore-5dm` (P3)
- [ ] **Duplicated code.** Remove the duplicated paths (`as_h5(NeuroVecSeq)`, `as.array.LatentNeuroVec`, `h5file` methods, `make_run_*` defers). `fmristore-skz` (P3)
- [ ] **Read performance.** Batch `linear_access` reads per cluster. `fmristore-0h2`

## Milestone 4: Readability and documentation hygiene

Ongoing; do opportunistically in the files each change touches.

- [ ] **Runnable examples.** Replace `\dontrun{}` examples with short tempfile-based examples that run. `fmristore-6z4` (P3)
- [ ] **One error and message style** (classed conditions for public APIs), enforced by lint. `fmristore-fub` (P3)
- [ ] **Roxygen.** Pin the roxygen2 version in CI and regenerate `man/` once. Several Rd files were edited by hand in this pass. `fmristore-5jx`
- [ ] **Orphans.** Remove the orphan roxygen header `#' Accessor Methods for H5ParcellatedMultiScan` in `R/cluster_experiment.R`. `fmristore-skz`

## Suggested order and size

| Order | Work | Size |
|-------|------|------|
| 0 | Green `main` (`fmristore-bfi`): fixture fix now, Windows and macOS triage | S-M |
| 1 | M0 (NEWS, CI pinning, lint/style gates, upstream report) | S-M |
| 2 | Decisions D1-D4 | - |
| 3 | M1 data-safety block (`cmi`, `qt9`, `71g`, `o4x`) | M |
| 4 | M1 broken-results block, one issue per PR | M-L |
| 5 | M1 documentation items, M4 examples for the same functions | S |
| 6 | M2 mapped into the 1.0 plan's phases and test suite | (1.0 schedule) |
| 7 | M3 per D1 | S (freeze) / L (full) |

Release a 0.x patch after steps 3-5.
