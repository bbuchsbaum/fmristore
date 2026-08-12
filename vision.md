# fmristore 1.0 vision

`fmristore` will become the transactional storage engine for `fmridataset`,
not a second neuroimaging data model.

Users will construct and analyze `fmri_frame`, `fmri_collection`, and
`fmri_study` objects. `fmristore` will make those objects durable, lazy,
versioned, verifiable, recoverable, and safely appendable.

## How this departs from the old package

| Old `fmristore` | Proposed `fmristore` 1.0 |
|---|---|
| Defines neuroimaging objects such as `H5NeuroVec`, latent vectors, experiments, groups, and parcellations | Stores semantic objects defined by `fmridataset` |
| Mixes domain semantics with HDF5 mechanics | Owns physical formats, transactions, integrity, recovery, and migration |
| Persistent objects may contain live HDF5 handles | Persistent descriptors are reconstructible and contain no live handles |
| Exposes a broad S4 class hierarchy | Exposes a small verb-based API plus `h5_array_source` |
| Embeds RDS-encoded manifests | Uses canonical, bounded, language-neutral manifests |
| Uses mutable root manifests without a safe concurrent-append protocol | Uses a SQLite transaction catalog with immutable HDF5 objects |
| Checks shape and metadata consistency | Supports layered manifest, object, and numerical-content verification |
| Uses paths and timestamps as part of source fingerprints | Separates semantic, content, storage, and locator identities |
| Acts partly as a general HDF5 convenience toolkit | Certifies two narrow storage formats |

## Two physical formats

The two formats serve different purposes:

- **FDS-H5 snapshot:** one immutable HDF5 file for publication, exchange,
  archival use, and fixed analysis inputs.
- **FDS directory store:** a SQLite catalog plus immutable HDF5 objects for
  incremental ingest, committed generations, and observation append.

The user-facing API remains small:

```r
create_store(frame, "analysis.fds.h5", mode = "snapshot")
frame <- open_store("analysis.fds.h5")

create_store(frame, "cohort.fds", mode = "appendable")
append_store("cohort.fds", next_subject)
verify_store("cohort.fds", level = "content")
```

The package will distinguish operations that have different safety and cost
profiles. Inspection reads metadata. Verification checks integrity. Recovery
handles interrupted transactions. Garbage collection removes unreferenced
immutable objects. None of these operations silently substitutes for another.

## Division of responsibility

`fmridataset` answers semantic questions: what an object means, how its axes
and identifiers align, and whether two objects may be combined. `fmristore`
answers physical questions: where arrays are stored, which generation refers
to them, whether their contents are intact, and how a write becomes visible
and durable.

This division is the central departure from the old design. The old package
asks, "What kind of neuroimaging object is this?" The new package asks, "How
can this already-defined scientific object be stored and recovered without
ambiguity or corruption?"

The change therefore requires a deliberate breaking reset. The legacy S4
hierarchy moves out of the 1.0 core. `h5_array_source`, bounded block I/O,
immutable data objects, inspection, migration, and failure-aware recovery form
the foundation of the new package.

The detailed scope, dependencies, phases, and release gates are in the
[fmristore 1.0 implementation plan](design/fmristore-1.0-plan.md).
