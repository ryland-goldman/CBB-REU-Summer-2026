# easygdf

Pure-Python library for reading and writing **GDF (General Datafile Format)** files —
the native binary I/O format of [General Particle Tracer (GPT)](http://www.pulsar.nl/gpt/).
A lightweight alternative to GPT's `GDF2A`/`ASCI2GDF` round-trip: load a `.gdf` straight
into Python/NumPy, modify, and save back. Source repository:
[easygdf on GitHub](https://github.com/electronsandstuff/easygdf).
BSD-3-Clause licensed; available on both PyPI (`pip install easygdf`) and conda-forge
(`conda install easygdf`).

The upstream README is the complete package reference (full signatures, parameter docs,
and the GDF block model); it lives at [source/README.md](source/README.md). The
`source/examples/` scripts are the practical worked reference.

---

## Overview

GDF files are organized into **blocks**, each a `dict` with three keys: `name` (ASCII
string), `value` (int / float / str / `None` / `bytes` / NumPy array — no complex types),
and `children` (a list of further blocks). `load` returns the whole file as a dict (header
fields + a `blocks` list); `save` accepts the same structure, so editing is a load →
mutate → save round-trip. For GPT's standard outputs, the convenience functions skip the
raw block tree and return clean per-screen / per-distribution dicts.

| Section | File |
|---------|------|
| Full package reference, quickstart, install, GDF block model | [source/README.md](source/README.md) |

---

## API reference (`import easygdf`)

| Function | Purpose |
|----------|---------|
| `load(f, ...)` | Read any GDF file into a dict (header fields + `blocks` list). |
| `save(f, blocks=None, ...)` | Write blocks (+ header) to a GDF file. Signature mirrors `load` output. |
| `load_screens_touts(f, ...)` | Parse a GPT output file into `{"screens": [...], "touts": [...]}` dicts. |
| `save_screens_touts(f, screens=None, touts=None, ...)` | Write a GPT-format output file; missing particle arrays auto-filled. |
| `load_initial_distribution(f, ...)` | Read a GPT initial-distribution file into a dict of arrays. |
| `save_initial_distribution(f, x=, y=, z=, GBx=..., ...)` | Write a GPT initial distribution; autofills required keys (zeros; sequential `ID`). |
| `is_gdf(f)` | Test whether a file/stream is a GDF file. |
| `get_example_screen_tout_filename()` | Path to a bundled example screens/touts `.gdf`. |
| `get_example_initial_distribution()` | Bundled example initial-distribution data. |
| `GDFError`, `GDFIOError` | Exception types. |
| `GDF_DOUBLE`, `GDF_INT32`, `GDF_ASCII`, … | GDF type-code constants. |

**Key field conventions** (see upstream README for the complete lists):
- Screen particle keys: `ID`, `x`, `y`, `z`, `Bx`, `By`, `Bz`, `t`, `m`, `q`, `nmacro`, `rmacro`, `rxy`, `G`
- Initial distribution: supply position (`x`,`y`,`z`) plus **either** momentum (`GBx`,`GBy`,`GBz`) **or** velocity (`Bx`,`By`,`Bz`) — not both.

---

## Examples (the practical reference)

| Topic | File |
|-------|------|
| Minimal raw-block read/write round-trip | [source/examples/minimal.py](source/examples/minimal.py) |
| Loading GPT screens & touts output | [source/examples/screens_touts.py](source/examples/screens_touts.py) |
| Writing a GPT initial particle distribution | [source/examples/initial_distribution.py](source/examples/initial_distribution.py) |
