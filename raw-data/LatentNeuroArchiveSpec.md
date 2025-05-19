Okay, this looks like a solid, self-contained specification document derived directly from our previous detailed discussion about the R package implementation. It correctly focuses on the *on-disk format* and *required semantics*, leaving implementation choices to the blueprint.

Here is the specification reformatted into Markdown for clarity and readability:

---

# Latent Neuro-Archive (LNA) Format Specification - Version 2.0

**Document Version:** 1.0 (Corresponds to R package blueprint v1.4)

This document provides the formal specification for the Latent Neuro-Archive (LNA) container format, version 2.0. It details the required on-disk structure within an HDF5 file and the semantics necessary for compliant reading and writing. Implementation details (specific classes, helper functions, etc.) are considered out of scope for this document and are covered by associated implementation blueprints (e.g., the `lna` R package v1.4 blueprint).

## 0. Scope & Terminology

| Term             | Meaning                                                                                                                                             |
| :--------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------- |
| **LNA file**     | A single HDF5 file, typically with extension `.lna.h5`, fully conformant to this document.                                                          |
| **Dataset**      | An HDF5 Dataset object (multidimensional array).                                                                                                    |
| **Group**        | An HDF5 Group object (container).                                                                                                                   |
| **Attribute**    | An HDF5 Attribute object (metadata attached to Groups or Datasets).                                                                                 |
| **Numeric Payload**| Any Dataset primarily storing large numerical arrays subject to transformation (especially quantisation), e.g., basis matrices, embeddings.         |
| **JSON Descriptor**| A Dataset of HDF5 type `H5T_C_S1` (variable-length, null-terminated, UTF-8 encoded string) containing a single JSON object describing one transform step. |
| **Transform**    | A logical forward step in the compression pipeline (e.g., quantisation, basis projection). Decoded by applying the inverse operation in reverse order. |
| **Core Transform**| One of the standard transforms defined by this spec: `quant`, `basis`, `embed`, `temporal`, `delta`.                                                 |
| **Required Transform**| A transform whose `type` string is listed in the root `Attribute` `required_transforms`. Readers MUST understand these to decode the file.        |
| **Voxel Indexing (N)**| Voxels within the `mask` are enumerated `0` to `N-1` based on their Fortran-order index (column-major) within the full `X x Y x Z` volume. `N = sum(mask)`. |

## 1. Versioning & Compatibility

*   The LNA format version is specified by the root `Attribute` `lna_spec = "MAJOR.MINOR"`.
    *   **Current Version: `"2.0"`**
*   A reader application supporting LNA version `S.T` MUST handle a file with `lna_spec` version `F.G` as follows:
    *   If `floor(F) > floor(S)`: The reader MUST refuse to open the file (incompatible major version).
    *   If `floor(F) == floor(S)` and `G > T`: The reader MAY open the file but SHOULD warn that the file uses newer optional features. It MUST ignore unknown optional transforms (§4).
    *   If `floor(F) < floor(S)`: The reader MAY open the file silently (assuming backward compatibility support).
*   Compatibility related to specific transform logic is handled via the `version` field within each transform's JSON Descriptor (§4.1).

## 2. Top-Level HDF5 Layout

The LNA file MUST contain the following HDF5 structure. Optional elements are noted.

```
/
├─ (Attributes)
│   ├─ lna_spec            : String Attribute ("2.0") - Required
│   ├─ creator             : String Attribute (Free text description of creating software) - Required
│   ├─ required_transforms : String Array Attribute (e.g., ["quant", "basis", "embed"]) - Required
│   └─ lna_checksum        : String Attribute (Optional, 64-char hex SHA-256 hash) - Optional (§3.5)
│
├─ header/              : Group - Required
│   └─ (Attributes)     : NIfTI-1 header fields stored as HDF5 Attributes - Required
│
├─ mask                 : Dataset {uint8, dims=[X,Y,Z]} (1=in-mask, 0=out) - Required
│
├─ spatial/             : Group - Required
│   ├─ partition            : Dataset {int32, dims=[N]} (Voxel index -> block_id; 0=global/none) - Required
│   ├─ block_table          : Dataset {int32, dims=[B, 7]} (Rows: x0,x1,y0,y1,z0,z1,k_b) - Required (if B > 0)
│   ├─ coeff_offset_block   : Dataset {int32, dims=[B]} (Column offset in embedding for block b) - Required (if B > 0)
│   └─ block_voxel_indices  : Dataset {int32, dims=[N]} (Optional: Voxel index -> local 0-based index within its block) - Optional
│
├─ basis/               : Group - Required (if basis transforms used)
│   ├─ global            : Dataset {numeric*, dims=[k_g, N] or [N, k_g]} (Global basis vectors) - Optional
│   └─ blocks/           : Group (Contains subgroups for block-specific bases) - Optional
│   │    └─ <id>/         : Group (Name is block_id, e.g., "1", "2", ...)
│   │        └─ matrix    : Dataset {numeric*, dims=[k_b, n_b] or [n_b, k_b]} (Basis vectors for block `b`)
│   └─ (Attributes)
│        └─ basis_reference : String Attribute (Optional: "sha256:<hash> url:<url>" for external basis) - Optional
│
├─ scans/               : Group - Required
│   └─ <run_id>/         : Group (One per run, name is run identifier, e.g., "run-01") - Required (at least one)
│        ├─ embedding      : Dataset {numeric*, dims=[T, k_total]} (Temporal coefficients) - Required
│        └─ (Attributes)
│             └─ metadata   : String Attribute (JSON string with run-specific info: TR, task, etc.) - Required
│
├─ plugins/             : Group (Storage for custom transform payloads) - Optional
│
└─ transforms/          : Group (Contains ordered JSON Descriptors) - Required
     ├─ 00_quant.json     : Dataset {H5T_STRING} (JSON Descriptor) - Required (typically)
     ├─ 01_basis.json     : Dataset {H5T_STRING} (JSON Descriptor) - Required (typically)
     ├─ 02_embed.json     : Dataset {H5T_STRING} (JSON Descriptor) - Required (typically)
     └─ ...               : Dataset {H5T_STRING} (Optional transforms) - Optional
```

*\*Note on `numeric*` types:* The exact HDF5 data type (e.g., `H5T_STD_I8LE`, `H5T_IEEE_F16LE`, `H5T_IEEE_F32LE`) for Numeric Payload datasets (`/basis/*`, `/scans/*/embedding`) is determined by the innermost `quant` transform applied to that payload. See §5 and §6.

## 3. Mandatory Semantics

### 3.1 Dimensions & Consistency

The following dimensional relationships MUST hold:

*   `sum(mask)` MUST equal `N`.
*   `N` MUST equal the number of columns (if `storage_order = "component_x_voxel"`) or rows (if `"voxel_x_component"`) in `/basis/global` if it exists.
*   All Datasets under `/spatial/` indexed by `N` must have `N` rows (e.g., `partition`, `block_voxel_indices`).
*   The `/header/` Attributes must reflect the original dimensions, e.g., `dim = [4, X, Y, Z, T_orig, 1, 1, 1]`.
*   `/scans/<run>/embedding` width `k_total` MUST equal `k_g` (number of global components) + sum of `k_b` for all blocks `b=1..B` specified in `/spatial/block_table`. The column ordering is global components first, followed by block components concatenated in ascending order of `block_id`.
*   For each block `b` defined in `/spatial/block_table`:
    *   `n_b` (the number of voxels within the block's extents *and* within the `mask`) MUST equal the number of columns (if `storage_order = "component_x_voxel"`) or rows (if `"voxel_x_component"`) in `/basis/blocks/<b>/matrix`.
    *   If `k_b = 0` for block `b`, the `/basis/blocks/<b>/` group MAY be omitted.

Readers MUST treat violations of these consistency rules as fatal format errors.

### 3.2 Nomenclature & ID Conventions

*   Block IDs (`block_id`) are integers starting from `1`. `block_id = 0` in `/spatial/partition` indicates the voxel uses only the global basis or is outside all defined blocks.
*   For a block with `block_id = b` (where `b >= 1`):
    *   Its properties are stored in row `b-1` (0-based index) of `/spatial/block_table`.
    *   Its specific basis matrix (if `k_b > 0`) is stored under the Group `/basis/blocks/<b>/` (where `<b>` is the string representation of the integer `b`).
    *   Its starting column offset within the block-specific part of the embedding matrix is given by `/spatial/coeff_offset_block[b-1]`.

### 3.3 JSON Descriptor Ordering

*   JSON Descriptor datasets under `/transforms/` MUST be named `NN_type.json`, where `NN` is a two-digit, zero-padded integer (`00` through `99`).
*   The sequence `NN` MUST be continuous, starting from `00`.
*   This sequence represents the *forward* order of applying transforms during compression (e.g., `00_quant` applied first, then `01_basis`).
*   Decoders MUST process transforms by inverting them in the *reverse* order of this sequence (e.g., invert `NN_...`, then `NN-1_...`, down to `00_...`).

### 3.4 Bit-Packing Rule (Sub-byte Integers)

*   When a `quant` transform specifies a logical integer type smaller than 8 bits (e.g., 4-bit or 6-bit), the data MUST be packed into bytes (`uint8`).
*   Packing MUST use little-endian order within the byte: the value for the first voxel index occupies the lowest-order bits (e.g., bits 0-3 for 4-bit data) of the first byte.

### 3.5 Checksum Attribute (`lna_checksum`)

*   This root `Attribute` is optional.
*   If present, its value MUST be a 64-character lowercase hexadecimal string representing the SHA-256 checksum of the *entire LNA file's byte stream*.
*   This checksum MUST be computed *after* the HDF5 file handle used for writing has been fully closed. (Storage may require reopening or external tooling).
*   Readers SHOULD verify this checksum if present. A mismatch SHOULD result in an error (suggested class: `"lna_error_checksum"`).

## 4. Transform System

### 4.1 JSON Descriptor Schema (Common Fields)

Every JSON Descriptor Dataset under `/transforms/` MUST contain a JSON object adhering to the following structure:

```json
{
  "type": "string",                // Required: Namespace.Type identifier (e.g., "quant", "acme.ae")
  "version": "string",             // Required: Version string for this transform's logic (e.g., "1.0")
  "params": { },                   // Required: Object containing transform-specific parameters. Structure defined by type/version. See Appendix B for core types.
  "datasets": [                    // Required: Array of objects referencing HDF5 datasets used/produced.
      {
        "path": "string",          // Required: Absolute HDF5 path to the dataset (e.g., "/basis/global").
        "role": "string"           // Required: Symbolic role (e.g., "basis_matrix", "coefficients", "lookup_table").
      }
  ],
  "inputs": ["string"],            // Required: Array of symbolic names (keys) expected in the decoder's intermediate state ('stash').
  "outputs": ["string"],           // Required: Array of symbolic names (keys) produced/replaced in the decoder's intermediate state ('stash').
  "metadata_updates": { },         // Optional: Object describing changes to canonical metadata fields (See Appendix A). Keys are field names, values are new effective values.
  "capabilities": {                // Optional: Object declaring support for partial decoding hints.
        "supports_spatial_subsetting": "boolean | 'block-only'", // Default: false
        "supports_temporal_subsetting": "boolean | 'contiguous-only'" // Default: false
  },
  "schema_uri": "string"           // Optional: URI pointing to a JSON Schema definition for the 'params' object (Recommended for non-core transforms).
}
```

*   **`inputs` / `outputs` Contract:** The decoding framework (e.g., the R package) MUST enforce that an inverse transform step consumes only the declared `inputs` from the intermediate state and produces exactly the declared `outputs`.
*   **`capabilities`:** Hints used by the decoder to potentially optimize reading when `roi_mask` or `time_idx` are provided. If a required transform lacks capability, the decoder MUST fall back to full decoding for that dimension.

### 4.2 Core Transforms

These transforms are defined as part of the LNA 2.0 standard. Readers SHOULD support them if encountered, and MUST support them if listed in `required_transforms`. See Appendix B for detailed parameter schemas.

| Type         | Forward Role Description                  | Inverse Requirement                                | Key Params (Illustrative)            |
| :----------- | :---------------------------------------- | :------------------------------------------------- | :----------------------------------- |
| `quant`      | Convert float Numeric Payloads -> int/f16 | De-quantise back to float                          | `mode`, `bits`, `scale`, `zero_point`  |
| `basis`      | Apply spatial basis projection            | Reconstruct voxel values or latent representation | `basis_type`, `storage_order`, `k_global` |
| `embed`      | Stack/select time courses into embedding  | Load correct `/scans/<run>/embedding` Dataset(s) | `run_selector`                     |
| `temporal`   | Apply temporal basis (DCT, B-spline...) | Multiply coefficients by inverse temporal basis | `kind`, `k`, `degree`, `knot_vector_path` |
| `delta`      | Delta or Run-Length Encode integers       | Undo differencing / RLE                            | `mode`, `order`, `reference`         |

### 4.3 Optional / Custom Transforms

Third-party or experimental transforms MUST use a namespaced `type` identifier (e.g., `"myorg.wavelet"`). They follow the same JSON Descriptor contract. Readers MAY skip unknown optional transforms if they are not listed in `required_transforms`.

## 5. Data Type Rules

Numeric Payloads and other datasets MUST adhere to these HDF5 type conventions based on their logical representation:

| Logical Type         | Recommended HDF5 Storage Type                    | Notes                                                                    |
| :------------------- | :----------------------------------------------- | :----------------------------------------------------------------------- |
| **float32**          | `H5T_IEEE_F32LE` / `H5T_IEEE_F32BE` (Native)     | Preferred final numeric type after decoding.                             |
| **float16**          | `H5T_IEEE_F16LE` / `H5T_IEEE_F16BE`              | Valid storage type. Reader typically converts to float32/64 upon load. |
| **int 4 / 6 / 8**    | Packed into `H5T_STD_U8LE` / `H5T_STD_I8LE`      | Use little-endian packing within bytes (§3.4). Signedness per `quant` params. |
| **int 16**           | `H5T_STD_I16LE` / `H5T_STD_I16BE`              | E.g., for `quant` mode `"linear_s16"`.                                 |
| **int 32 / uint32**  | Native signed/unsigned 32-bit integer          | E.g., for `partition`, `block_table`, indices.                         |
| **uint8**            | `H5T_STD_U8LE`                                   | E.g., for `mask`.                                                        |

*   **Endianness:** Writers SHOULD use native endian types (e.g., `H5T_NATIVE_FLOAT`). HDF5 handles translation between file (typically LE) and machine endianness for standard types. Sub-byte packing endianness is defined in §3.4.
*   **Determination:** The exact stored HDF5 type for a Numeric Payload is determined by inspecting the Dataset's properties (`H5Tget_class()`, `H5Tget_size()`, `H5Tget_sign()`, `H5Tget_order()`). The `quant` transform parameters define the *logical* interpretation and conversion.

## 6. Validation Requirements (Non-Exhaustive)

A conformant LNA reader MUST perform validation checks including, but not limited to:

*   Presence and readability of all mandatory HDF5 Groups, Datasets, and Attributes (§2).
*   `lna_spec` compatibility (§1).
*   Transform chain numbering continuity and format (`NN_type.json`) (§3.3).
*   Existence of all HDF5 paths listed in the `datasets` array of every *required* transform's JSON Descriptor.
*   Basic HDF5 data type class compatibility for paths listed in `datasets` (e.g., `quant` expects numeric, not string).
*   Dimensional congruence rules specified in §3.1.
*   Verification of `lna_checksum` if present and validation requested (§3.5). Error class: `"lna_error_checksum"`.
*   Availability of an implementation for the inverse operation of every transform listed in the `required_transforms` attribute. Error class: `"lna_error_unknown_required_transform"`.
*   Validity of `params` object against the known JSON schema for core transform types (Appendix B).

Failure of any MUST-level requirement during validation implies the file is non-conformant or cannot be decoded by the current reader.

## 7. Appendices

### Appendix A: Canonical Metadata Keys

List of standard keys recognized by the `metadata_updates` field in JSON Descriptors. Keys typically align with NIfTI-1 / BIDS standards. (Units apply to the *value* provided in `metadata_updates`).

*   `dim`: (int vector) Array dimensions.
*   `pixdim`: (numeric vector) NIfTI pixdim array.
*   `TR`: (numeric, **seconds**) Effective repetition time.
*   `voxel_spacing`: (numeric vector, units from original `pixdim`/`xyzt_units`) `c(dx, dy, dz)`.
*   `origin`: (numeric vector, units from original `pixdim`/`xyzt_units`) Voxel coordinates of origin.
*   `axes_orientation`: (string) e.g., "RAS", "LPI".
*   `time_units`: (string) e.g., `"seconds"`, `"msec"`, `"hz"`.
*   `spatial_units`: (string) e.g., `"mm"`, `"micron"`.
*   `qfac`: (numeric, +/- 1) From `pixdim[0]`.
*   *(This list can be extended in minor revisions of the spec).*

### Appendix B: Core Transform Parameter JSON Schemas

The detailed JSON Schema definitions for the `params` object of the core transforms (`quant`, `basis`, `embed`, `temporal`, `delta`) are maintained alongside the reference implementation (e.g., in the `inst/schemas/` directory of the `lna` R package). These schemas define required fields, data types, allowed values, and default values for parameters.

### Appendix C: Octree Spatial Transform Example (Illustrative)

This shows how a non-core, block-based spatial transform might be described. It implies a specific internal structure (e.g., under `/plugins/` or within `/basis/`) not mandated by the core spec.

```json
{
  "type": "myorg.octree_basis",
  "version": "0.1",
  "params": {
      "depth": 3,                     // Max octree depth
      "min_components_per_node": 10,  // Minimum k per node
      "storage_order": "component_x_voxel"
  },
  "datasets": [ // Example paths - actual paths defined by transform logic
      { "path": "/plugins/myorg_octree/node_basis_matrices", "role": "basis_collection" },
      { "path": "/plugins/myorg_octree/voxel_to_node_map",   "role": "partition_map" },
      { "path": "/plugins/myorg_octree/node_offset_table",   "role": "offset_table" }
  ],
  "inputs":  ["embedding"],          // Consumes coefficients
  "outputs": ["latent"],             // Produces reconstructed latent data
  "capabilities": {
      "supports_spatial_subsetting": true, // Can prune tree traversal
      "supports_temporal_subsetting": true
  },
  "schema_uri": "https://myorg.com/schemas/lna/octree_basis_v0.1.schema.json"
}
```
*Reader Handling:* A reader supporting `myorg.octree_basis` would use the specified `datasets` (e.g., `partition_map` mapping voxels to octree node IDs, `basis_collection` holding matrices, `offset_table` for embedding column lookups) to perform the inverse basis multiplication, potentially pruning the operation based on spatial subset hints.

## 8. Conformance Statement

*   An **LNA file** conforms to LNA 2.0 if it meets all structural and semantic requirements outlined in Sections 1 through 6 of this document.
*   An **LNA reader** conforms to LNA 2.0 if it:
    1.  Correctly validates files against the requirements in Sections 1-6.
    2.  Accurately implements the inverse operation for all core transforms (`quant`, `basis`, `embed`, `temporal`, `delta`) as defined by this specification and their respective parameter schemas (Appendix B), when they are listed as required or encountered.
    3.  Can correctly skip unknown optional transforms (those not listed in `required_transforms`) without compromising the decoding correctness of subsequent required transforms.

---
**END OF SPECIFICATION**