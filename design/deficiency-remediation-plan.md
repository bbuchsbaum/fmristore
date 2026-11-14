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
- **Legacy removal, step 2 (done):** delete the remaining parcellated and
  labeled classes. neurotabs and fmrigds are consumers and move to frames on
  their own schedule.
- **1.0 engine requirements:** frame-store durability and integrity findings
  become requirements and storage-law tests of the new engine.

## Decisions

| ID | Decision | Consequence |
|----|----------|-------------|
| D1 | Remove legacy code now, as long as fmridataset is unaffected (it uses only the frame store). fmristore is self-contained; neurotabs and fmrigds are consumers and adapt to it. | Step 1 deleted everything no package calls. Step 2 (0.2.0) deleted the rest with no converters, since no real users hold legacy files. Consumers can pin to `4d7adb7` until they migrate. |
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
- [ ] **Pin dependencies and test the installed package in CI.** Pin `Remotes:` to commits or tags. `fmristore-5jx` (P2)
- [ ] **Lint gate.** Add `any::cyclocomp` to `lint.yaml`, drop the two linters per D2, and lint only changed files. Exit: the lint job reports and passes. `fmristore-ky0` (P2)
- [ ] **Style gate.** Check style on changed files only. Exit: the style job passes. `fmristore-ky0`
- [x] ~~**Upstream `concat` bug.**~~ Out of scope: fmristore no longer handles latent or neuroim2 vectors. `fmristore-qsh` (closed)

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

## Milestone 2: Legacy removal, step 2 (done)

- [x] **Deleted** the parcellated and labeled classes, `read_dataset()`,
  `detect_h5_type()`, the example helpers, constants and memory settings
  (about 5,000 lines). Imports are now fmridataset and hdf5r. `fmristore-w1f`
- [x] Converters were written and then dropped: no real users hold legacy
  files. `fmristore-qv9`
- [x] Frozen-class bugs closed as moot: `fmristore-o4x`, `fmristore-e2u`,
  `fmristore-qt9`. The `volume_space` to neuroim2 helper is not fmristore's
  job (`fmristore-amb`, closed).

Consumer work, in those repositories (tracked here for reference):

- [ ] **neurotabs:** store parcel data as frames. `fmristore-dtd`
- [ ] **fmrigds:** regenerate fixtures with `write_frame_h5()`, detect with
  `frame_store_info()`. `fmristore-au7`

How their needs map to frames:

| Need | Frame store |
|---|---|
| Time x parcel summary per run | `run_id` observation column plus `filter_observations()` |
| One row by index | `source_read(observations = i)` |
| Parcel IDs, voxel counts, mask | `parcel_space`, `parcel_membership()`, parent `volume_space` |
| Named maps under a mask | observations = labels, features = mask voxels |
| File type detection | `frame_store_info()` |

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
| 2 | M2 (done); neurotabs and fmrigds migrate in their repos | - |
| 3 | M3 inside the 1.0 phases | (1.0 schedule) |

M4 runs alongside.
