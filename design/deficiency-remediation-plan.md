# Deficiency remediation plan

Status: decisions recorded 2026-11-13; removal step 1 done 2026-11-14. Tracks the deficiencies found while
fixing the test suite, reviewing `R/`, and rewriting the vignettes. Every item
links to a beads issue in `.beads/issues.jsonl`.

## How this relates to the 1.0 plan

`design/fmristore-1.0-plan.md` replaces the frame store with a new engine
(Phases 1-5) and removes the legacy surface (`H5Neuro*`, `H5Parcellated*`,
`LabeledVolumeSet`, latent I/O) in Phase 6. The legacy classes have no real
users outside our own packages, so we do not maintain them (D1). Work splits
into:

- **Green CI and an honest changelog** (Milestone 0).
- **Legacy removal, step 1 (done):** delete everything no package calls.
- **Legacy removal, step 2:** move neurotabs and fmrigds to frames, then
  delete the remaining parcellated and labeled classes.
- **1.0 engine requirements:** frame-store durability and integrity findings
  become requirements and storage-law tests of the new engine.

## Decisions

| ID | Decision | Consequence |
|----|----------|-------------|
| D1 | Remove legacy code now, in two steps, as long as fmridataset is unaffected (it uses only the frame store). | Step 1 deleted everything no package calls. Step 2 deletes the parcellated and labeled classes once neurotabs and fmrigds move to frames. |
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
  - [x] Ubuntu: fix the malformed fixture at `test-h5_corruption_handling.R:122`. `matrix(rnorm(8 * 5), nrow = 5, ncol = 10)` and the `7 * 5` line under-fill the matrix, which R 4.6 makes an error. The test runs only with `NOT_CRAN=true`.
  - [ ] Windows (11 failures): triage from the check artifact. The visible one is `test-vec-regressions.R:278` (`H5NeuroVec(out)`), likely an open handle blocking a copy or reopen. Related: `fmristore-57h`, `fmristore-71g`.
  - [ ] macOS: triage from the check artifact.
  - [x] Ubuntu and Windows: the upstream `concat` failure. The test went with latent I/O in step 1; the bug still needs reporting upstream (`fmristore-qsh`).
  - [x] Run local and agent checks with `NOT_CRAN=true` and `_R_CHECK_MATRIX_DATA_=TRUE`, as CI does (recorded in CLAUDE.md).
  - [ ] Confirm on the next hosted run; most Windows failures were in code removed by step 1.
- [x] **Changelog.** NEWS.md records the merged behaviour changes and the step 1 removals. `fmristore-2v0`
- [ ] **Pin dependencies and test the installed package in CI.** Pin `Remotes:` to commits or tags. Add a test convention for comparing neuroim2 objects via `@.Data`. `fmristore-5jx` (P2)
- [ ] **Lint gate.** Add `any::cyclocomp` to `lint.yaml`, drop the two linters per D2, and lint only changed files. Exit: the lint job reports and passes. `fmristore-ky0` (P2)
- [ ] **Style gate.** Check style on changed files only. Exit: the style job passes. `fmristore-ky0`
- [ ] **Upstream `concat` bug.** File it on fmrilatent (or neuroim2's `.concat4D`) with the reproduction and a test there. The fmristore test was deleted with latent I/O. `fmristore-qsh` (P3)

Exit evidence: a green hosted R CMD check matrix on one pinned dependency set.

## Milestone 1: Legacy removal, step 1 (done)

Deleted everything no package calls: `H5NeuroVol`, `H5NeuroVec`,
`H5NeuroVecSeq`, `to_nih5_vec()`, `neurovecseq_to_h5()`, latent I/O and the
fmrilatent re-exports, `as_h5()`, `write_dataset()`, `write_cluster_result()`,
`write_parcellated_scan_h5()`, `write_clustered_experiment_h5()` and
`make_run_*()`, with their tests, vignettes, specs and fixtures. fmristore no
longer imports fmrilatent. Exports went from 55 to 37.

- [x] fmridataset checked: it calls only `write_frame_h5()`, `open_frame_h5()` and `h5_array_source()`.
- [x] Issues in removed code closed: `cmi`, `bpi`, `crp`, `x7o`, `57h`, `71g`.
- [x] Docs regenerated with roxygen2 7.3.3 (the pinned `RoxygenNote`).
- [x] One vignette on storing frames replaces the five legacy vignettes.

## Milestone 2: Legacy removal, step 2

Still in use, and frozen until removed:

- **neurotabs:** `read_dataset()`, `H5ParcellatedMultiScan` and `H5ParcellatedArray` (`@runs`, `@clusters`, `@mask`, `as.matrix()`), and `write_parcellated_experiment_h5()` in tests.
- **fmrigds:** `detect_h5_type()`, and in tests `write_labeled_vec()` and `create_minimal_h5_for_H5ParcellatedMultiScan()`.

A prototype confirmed the frame store covers both packages' needs:

| Need | Frame store |
|---|---|
| Time x parcel summary per run | `run_id` observation column plus `filter_observations()` |
| One row by index | `source_read(observations = i)` |
| Parcel IDs, voxel counts, mask | `parcel_space`, `parcel_membership()`, parent `volume_space` |
| Named maps under a mask | observations = labels, features = mask voxels |
| File type detection | `frame_store_info()` |

Checklist:

- [ ] **Converters** from parcellated and labeled files to frames, written while the legacy readers still exist. `fmristore-qv9` (P1)
- [ ] **`volume_space` to `neuroim2`** helper in fmridataset, for masks, centroids and space names. `fmristore-amb` (P2)
- [ ] **neurotabs:** `fmristore-parcel` backend and parcel ingest read frames. `fmristore-dtd` (P2)
- [ ] **fmrigds:** adapter and fixtures use frames; regenerate `labeled_volume_small.h5`. `fmristore-au7` (P2)
- [ ] **Delete** the remaining classes, readers and writers, and the `neuroim2`, `Matrix` and `assertthat` imports. `fmristore-w1f` (P2)
- [ ] Until then, fix a frozen-class bug only if it breaks CI or blocks neurotabs or fmrigds: `fmristore-o4x` (label 0), `fmristore-e2u` (`LabeledVolumeSet` space), `fmristore-qt9` (`write_labeled_vec()` overwrite).

Caveat: reading the voxels of one cluster from a voxel-level frame selects
scattered columns. If neurotabs or fmrigds need that path, order features by
cluster or choose a matching chunk layout.

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
| 0 | Confirm green `main` after step 1 (`fmristore-bfi`) | S |
| 1 | Rest of M0: CI pinning, lint and style gates, upstream `concat` report | S-M |
| 2 | M2 converters and the `volume_space` helper | M |
| 3 | M2 neurotabs and fmrigds migrations (in those repos) | M |
| 4 | M2 deletion | S |
| 5 | M3 inside the 1.0 phases | (1.0 schedule) |

M4 runs alongside.
