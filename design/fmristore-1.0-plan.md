# fmristore 1.0 implementation plan

- Status: Proposed
- Updated: 2026-08-12
- Planning baseline: `7874717`

## Purpose

`fmristore` 1.0 will persist canonical `fmridataset` objects. It will own
physical formats, array bindings, transactions, integrity checks, recovery,
and migration. It will not define a second neuroimaging object model.

The package will support two physical formats:

1. An immutable HDF5 snapshot for exchange, archival use, and fixed analysis
   inputs.
2. An appendable directory store with a SQLite catalog and immutable HDF5 data
   objects.

Development can begin before `fmridataset` stabilizes, but the logical schema,
canonical manifest encoding, semantic digest, and append semantics cannot be
frozen until the upstream contract is released. Work done before then must sit
behind a narrow internal adapter and use synthetic conformance fixtures. It
must not copy semantic binding rules into `fmristore`.

## 1.0 outcome

An ordinary user should work with `fmri_frame`, `fmri_collection`, or
`fmri_study` objects and call store verbs. The minimum 1.0 API is:

```r
create_store(x, path, mode = c("snapshot", "appendable"), ...)

open_store(
  path,
  generation = "latest",
  verify = c("manifest", "none", "objects", "content")
)

append_store(path, x, expected_generation = NULL, ...)
replace_store(x, path, ...)

inspect_store(path)
verify_store(path, level = c("manifest", "objects", "content"))
recover_store(path, apply = FALSE)
gc_store(path, apply = FALSE)
migrate_store(path, destination, replace = FALSE)
snapshot_store(path, destination, generation = "latest")
```

`h5_array_source` remains public as the reconstructible HDF5 binding used by
the `fmridataset` array-source protocol. A read-planning helper may also remain
public for diagnostics. Live HDF5 and SQLite handles remain internal and must
never be serialized into a frame.

The required 1.0 object scope is an `fmri_frame` snapshot plus observation
append to an `fmri_frame` directory store. Collection and study snapshots may
join 1.0 only if their upstream manifests pass the same conformance gates.
Collection or study mutation is not required for 1.0.

Inspection, verification, recovery, and garbage collection return small
serializable S3 reports such as `store_info`, `verify_report`, and
`recovery_report`. Structured conditions derive from `fmristore_error` and
carry the operation, store path, generation, object ID, and expected and actual
values when those fields apply. Callers must not have to parse printed text to
handle a conflict or integrity failure.

### Required guarantees

- Opening a store reads manifests and physical metadata, not numerical array
  values.
- A snapshot is immutable after publication.
- A directory-store reader observes one committed generation until reopened.
- A directory-store append commits one complete generation or no generation.
- An optional expected generation detects concurrent modification.
- Numerical content has a stable identity that does not include its locator.
- Copying a store changes its locator but not its semantic or content digests.
- Writers, readers, and verifiers obey documented memory budgets.
- Inspection, verification, recovery, and garbage collection are separate
  operations with separate reports.
- Durability claims name the operating systems and filesystem profiles on
  which they were tested.

## Package boundary

`fmridataset` owns:

- frame, collection, and study semantics;
- axes, IDs, feature spaces, maps, entities, relations, and tables;
- source-free logical manifests and semantic digests;
- schema compatibility and observation-bind planning;
- reconstruction from a manifest and physical bindings.

`fmristore` owns:

- snapshot and directory format identifiers;
- HDF5 dataset layout, chunks, filters, and byte order;
- physical dtype conversion and validity encoding;
- immutable object identifiers and locators;
- content and storage digests;
- atomic publication and durable commit;
- catalog transactions, generations, and key indexes;
- verification, recovery, garbage collection, and migration.

`fmristore` must call supported `fmridataset` functions at this boundary. It
must not mutate frame fields or infer semantic compatibility from list shape.

The expected core dependencies are `fmridataset`, `hdf5r`, `DBI`, and
`RSQLite`, plus the selected streaming-hash and canonical-manifest
implementation. `Matrix` remains only if a certified sparse physical codec
needs it. The core engine should not import `neuroim2` or `fmrilatent`.
Dependencies are selected and added only after their design spikes establish
the required behavior.

## Versioned contracts

The implementation must record three independent identities:

1. The logical schema and version, owned by `fmridataset`.
2. The physical format and version, owned by `fmristore`.
3. The writer implementation and version.

The proposed physical identifiers are:

```text
org.fmristore.fds-h5/1
org.fmristore.fds-dir/1
```

The final logical identifier is an upstream decision. Existing
R-serialization layouts must be treated as provisional read-only inputs even
when their current attributes contain `v1`.

## Physical formats

### FDS-H5 snapshot

The snapshot is one immutable HDF5 file. It contains:

- format, logical-schema, object-type, and commit metadata;
- one bounded canonical manifest;
- generated physical dataset names rather than semantic names;
- declared physical dtype, chunk grid, filters, and physical-axis order;
- content and storage digests for every numerical object;
- no R-serialized object graph and no absolute locator in object identity.

`open_store()` recognizes this format from its declared magic and identifiers,
not from its filename extension.

`create_store(..., mode = "snapshot")` errors when the destination exists.
`replace_store()` is the only public operation that replaces an existing
snapshot.

### FDS directory store

The appendable layout is:

```text
cohort.fds/
|-- FORMAT
|-- catalog.sqlite
|-- objects/
`-- staging/
```

HDF5 files under `objects/` are immutable. The SQLite catalog stores:

- store identity and format version;
- committed generations and parent generations;
- canonical manifest bytes and hashes;
- object declarations and relative locations;
- array bindings;
- observation and entity key indexes;
- transaction records and schema migrations.

`open_store()` recognizes a directory store from `FORMAT` and validates the
catalog identity before reading generation metadata.

The catalog is a transactional projection of the canonical manifest. It must
not become a competing semantic model.

The initial concurrency contract is one serialized commit at a time with any
number of fixed-generation readers. Local filesystems are the first certified
profile. Network filesystems require one externally serialized writer until a
dedicated conformance suite proves stronger behavior.

## Identity and integrity

The format distinguishes four concepts:

- **Semantic digest:** logical meaning, supplied by `fmridataset`.
- **Content digest:** canonical logical numerical values and dtype.
- **Storage digest:** physical declaration, including chunks, filters,
  physical-axis order, and content digest.
- **Locator revision:** the current relative path or URI. It is not identity.

The content digest must not depend on writer block size or HDF5 chunk layout.
The format specification must define a canonical value order or a fixed
logical-tile Merkle construction. The writer computes this digest while it
streams converted values; it must not perform a second full read solely for
hashing.

Verification levels are cumulative:

```text
manifest  decode bounded canonical bytes and verify schema and manifest hash
objects   verify referenced objects, datasets, dtypes, chunks, and filters
content   stream logical values and verify content digests
```

`inspect_store()` reports commit state and the last verified level separately.
It must not label a store complete merely because a root marker is present.

## Dtype and missing-value contract

Every binding records:

```text
logical_dtype
storage_dtype
realized_dtype
```

The first certified profile should support only types that have an exact,
tested R conversion:

```text
uint8  int8  uint16  int16  uint32  int32  float32  float64
```

Logical values require an explicit value encoding and validity encoding so
that missing values are not confused with ordinary NaNs or integer sentinels.
Logical arrays may use a byte value array plus a validity bitmap.

Defer `float16`, `bfloat16`, complex values, and full-range `int64` and
`uint64` until their conversion, missingness, and realized R representation
have written round-trip laws.

Read planning reports:

- physical bytes touched;
- realized R bytes returned;
- estimated peak working bytes.

## Transaction protocol

An append uses this sequence:

1. Create a UUID transaction record containing the base generation and writer
   identity.
2. Stream each new HDF5 object to the transaction staging directory.
3. Convert values and compute content and storage digests while writing.
4. Flush and close HDF5, synchronize the file, publish it under a unique
   immutable object name, and synchronize its containing directory.
5. Begin a short SQLite transaction.
6. Recheck the expected generation.
7. Insert object references, observation and entity keys, and array bindings.
8. Let database uniqueness constraints reject key collisions.
9. Insert the new generation and canonical manifest hash.
10. Mark the transaction committed and commit SQLite.
11. Remove staging metadata.

A crash before the catalog commit leaves no visible generation. Published but
unreferenced objects remain available for explicit garbage collection. A crash
after the catalog commit leaves a complete generation whose objects were made
durable first.

## Recovery and garbage collection

Each transaction record contains:

- transaction ID;
- writer host, PID, and process-start token;
- base generation;
- start time and lease state;
- staged and published object IDs;
- intended commit and transaction state.

Recovery distinguishes active, stale, committed, and abandoned transactions.
It removes only artifacts owned by a transaction proven stale or abandoned.
Garbage collection handles durable unreferenced objects. Verification handles
missing or corrupt committed objects.

Every deletion target is resolved relative to the store root. The
implementation rejects absolute paths, parent traversal, symlink escapes, and
unexpected object names.

## Work plan

### Phase 0: Freeze the legacy baseline

Can start now.

Deliverables:

- Tag and preserve the final legacy implementation and fixture set.
- Stop describing the RDS-backed layouts as final FDS v1 formats.
- Add read-only format recognition for the current snapshot, shard, and
  `0.1-provisional` layouts.
- Inventory internal and public consumers of every exported legacy class and
  function.
- Decide which legacy layouts receive migrators and which receive documented
  refusal.
- Pin CI to immutable dependency revisions instead of default development
  branches. Release builds later use tagged minimum versions.
- Record decisions for the package boundary, manifest encoding, SQLite
  journal and durability profile, and legacy compatibility policy.
- Restore green package check, style, lint, and focused codec workflows.

Exit evidence:

- A clean hosted R CMD check matrix from one exact dependency set.
- A consumer migration table with an owner and replacement for every retained
  use.
- Golden legacy fixtures stored independently of generated tests.

### Phase 1: Isolate the new engine

Can start now.

Deliverables:

- Introduce the compact store API and structured condition hierarchy.
- Move new code behind storage-specific modules for sources, snapshots,
  catalogs, transactions, integrity, recovery, and migration.
- Harden `h5_array_source` so it validates physical datatype, byte order,
  chunks, and filters rather than trusting attributes.
- Use store-relative object references where possible.
- Make resource ownership mechanical with scoped HDF5 and SQLite helpers.
- Define one internal semantic-adapter interface and test it with synthetic
  manifests. Do not publish or freeze that adapter yet.

Exit evidence:

- No new engine object contains a live persistent handle.
- Repeated failed opens, missing datasets, double close, and read-after-close
  leave no resources open.
- The engine can be tested without loading the legacy S4 hierarchy.

### Phase 2: Complete format and durability prototypes

Can start now, but remains provisional until Phase 3.

Deliverables:

- Draft snapshot, directory, canonical-manifest, dtype, and integrity
  specifications under `inst/spec/`.
- Compare canonical CBOR with a strictly typed canonical JSON profile using R
  and Python proof readers.
- Implement bounded decoder prototypes with byte, depth, field-count, and
  string-length limits.
- Implement chunk-independent streaming content hashes.
- Add native helpers for file synchronization, directory synchronization, and
  atomic replacement on supported platforms.
- Prototype the SQLite schema and transaction state machine with synthetic
  semantic plans.

Exit evidence:

- R and Python produce identical bytes and hashes for synthetic golden
  manifests.
- Process termination at each publication boundary yields only an old or new
  snapshot.
- The catalog state machine passes crash and generation-conflict tests without
  depending on unstable frame internals.

### Phase 3: Bind the certified snapshot to fmridataset

Wait for a tagged stable `fmridataset` contract.

Required upstream functions or equivalents:

- source-free frame and study manifests;
- canonical semantic digests;
- physical binding extraction;
- supported reconstruction from a manifest and bindings;
- binding validation without direct field access.

Deliverables:

- Pin the minimum compatible `fmridataset` version.
- Finalize the logical-to-wire manifest mapping and format identifiers.
- Implement snapshot creation, replacement, opening, inspection, and layered
  verification.
- Make opening metadata-only and relocation-safe.
- Support `fmri_frame` snapshots as the required object type. Add study and
  collection snapshots only if their upstream manifests meet the same gate.

Exit evidence:

- Semantic and numerical round trips pass from fresh installed packages.
- R and Python read the same golden snapshots.
- Corrupt manifests, values, datatype attributes, chunks, and missing objects
  fail at the documented verification level.

### Phase 4: Bind transactional append to fmridataset

Wait for stable upstream append planning.

Required upstream functions or equivalents:

- schema comparison;
- keyed observation-bind planning;
- entity, relation, table, and provenance merge planning;
- composition of lazy row partitions.

Deliverables:

- Implement `create_store(..., mode = "appendable")`.
- Implement expected-generation append and fixed-generation open.
- Add indexed observation and entity keys with uniqueness constraints.
- Execute the upstream semantic plan without recreating its rules.
- Keep append metadata work proportional to new rows and indexed lookups.

Exit evidence:

- Two-writer and eight-writer barrier tests show no lost updates.
- Duplicate observation or conflicting entity IDs commit at most once.
- Readers see one internally consistent generation.
- Existing committed HDF5 objects remain byte-identical after append.

### Phase 5: Finish recovery and migration

Starts after the directory transaction model is stable.

Deliverables:

- Implement transaction-aware recovery and explicit garbage collection.
- Migrate the current RDS snapshot and shard layouts.
- Migrate selected legacy neuroimaging layouts through supported semantic
  objects rather than recreating the old public classes.
- Make replacement and in-place migration failure-safe.
- Add `snapshot_store()` for exporting a committed directory generation.

Exit evidence:

- Every released or provisional fixture opens, migrates, or returns a
  documented structured refusal.
- Active writer artifacts are never classified as removable.
- Failed migrations preserve the source and any prior destination generation.

### Phase 6: Remove the old package surface and migrate consumers

Starts when the replacement snapshot API is certified.

Known consumers to address include:

- the legacy HDF5 backend in `fmridataset`;
- parcellated ingestion in `neurotabs`;
- labeled, latent, and parcellated adapters and fixtures in `fmrigds`;
- the Python ports that mirror the old class hierarchy.

Deliverables:

- Move necessary legacy decoders into a read-only companion or internal
  migration module.
- Remove `H5Neuro*`, `H5Parcellated*`, `LabeledVolumeSet`, latent object I/O,
  generic HDF5 convenience exports, and format guessing from the main package.
- Remove hard `neuroim2` and `fmrilatent` imports from the core engine.
- Update DESCRIPTION, README, NEWS, vignettes, and the pkgdown reference around
  the store API and its guarantees.
- Convert one Python implementation into the cross-language conformance
  reader; freeze the other legacy port rather than duplicating new work.

Exit evidence:

- No supported consumer calls a removed class or generic.
- The installed namespace exposes only the store API, source protocol methods,
  reports, and conditions intended for 1.0.
- Core package checks do not require downstream development packages.

### Phase 7: Certify and release 1.0

Deliverables:

- Run core package checks on R devel, release, and oldrel across Linux, macOS,
  and Windows.
- Run snapshot conformance, dtype, relocation, resource, and migration jobs.
- Run real process-death and competing-writer jobs.
- Run downstream `fmridataset` and `fmrigds` integration separately from the
  core package gate.
- Certify local filesystem behavior and document unverified filesystem
  profiles.
- Publish specifications and golden stores for every supported format.

Release gate:

- Zero R CMD check errors, warnings, or unexplained notes.
- No default-branch dependency resolution in release workflows.
- All required storage laws pass on every supported platform.
- Documentation states the exact atomicity, durability, concurrency, and
  verification guarantees.
- Every provisional physical layout has a tested migration route.
- The physical format is frozen only after the release candidate passes the
  full conformance matrix.

## Storage-law test suite

Tests are organized around guarantees rather than source files alone.

### Atomicity

Terminate writes during each array, after HDF5 flush and close, before and
after object publication, during the catalog transaction, and after catalog
commit. Reopening must yield the complete old generation or the complete new
generation.

### Process death

Run writers in separate processes and terminate them with uncaught errors,
normal immediate exit, SIGTERM, and SIGKILL where supported. An injected R
error is not sufficient because `on.exit()` still runs.

### Concurrency

Use process barriers for simultaneous writers. Test generation conflicts,
duplicate keys, unique object names, and fixed-generation readers.

### Integrity

Corrupt manifest bytes, values, datatype attributes, chunks, object files,
generation parents, and referenced paths. Each mutation must fail at its
documented verification level.

### Dtypes and missingness

Test finite limits, zero and negative zero, missing values, ordinary NaNs,
infinities, unsigned limits, byte order, rejected overflow, and declared
round-trip tolerances.

### Relocation and resources

Copy stores to new paths and verify stable identities and lazy reads. Repeated
failure must leave no HDF5 or SQLite handles open, including immediate file
deletion tests on Windows.

### Memory and scale

Measure peak memory for writers and verifiers against the documented budget
plus bounded codec overhead. Scale fixtures must cover many assays, metadata
blocks, immutable objects, and at least hundreds of thousands of observation
keys without making open or append scan existing numerical values.

## Explicit non-goals for 1.0

- A general HDF5 object mapper or arbitrary HDF5 reader.
- Mutable committed HDF5 datasets.
- HDF5 SWMR as the transaction mechanism.
- Cloud object stores, Zarr, or remote catalogs.
- Transparent multiwriter guarantees on untested network filesystems.
- Plugin-only compression as the required default.
- Automatic repair of corrupt committed values without a redundant copy.
- Neuroimaging image, latent, parcellation, experiment, or group classes.
- Statistical analysis APIs.

## Immediate work while fmridataset stabilizes

The first implementation sequence is:

1. Freeze and label the legacy formats, pin dependencies, and restore a clean
   baseline.
2. Isolate the new engine, condition hierarchy, resource ownership, and
   hardened `h5_array_source`.
3. Draft the physical specifications and complete the canonical-encoding and
   durability spikes.
4. Build the SQLite catalog and transaction state machine against synthetic
   semantic fixtures.
5. Stop at the semantic adapter. Do not freeze manifests, digests, append
   compatibility, or public 1.0 format identifiers until a tagged
   `fmridataset` contract is available.

This sequence produces useful, testable storage machinery without embedding a
temporary version of `fmridataset` semantics into the physical engine.
