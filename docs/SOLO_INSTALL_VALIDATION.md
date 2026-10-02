# Solo installation validation

Checked on 2026-10-01 on macOS ARM64 with Apple Clang 21.0.0 and CMake 4.3.4.

| Check | Result |
|---|---|
| Make build and installation | Passed |
| CMake Release build and installation | Passed |
| `bash tests/run_tests.sh` | 47 passed, 0 failed |
| Installed-tool smoke test, Make and CMake | Passed |
| CMake-built 64-bit index, read back by installed inspect | Passed |
| Real EmptyDrops helper with NumPy 2.4.6 | Completed for Gene and GeneFull |
| Missing NumPy | Visible diagnostic; knee matrices preserved |

The installed-tool smoke test uses an installation, index and output directory
containing spaces and apostrophes, and runs from an unrelated working directory.
It checks the wrappers, gene-model extractor, both gene matrices, allele output,
helper lookup and failure fallback. Helper spies isolate lookup from statistical
behaviour. The separate real NumPy check used a small synthetic library with no
eligible rescue candidates; it verifies execution, not cell-calling accuracy.

The regression suite checks known REF/ALT counts, gene counting, thread
reproducibility, index construction and other alignment behaviour. These tests
do not establish equivalence to CellRanger on real human samples.

The Rust builder remains a separate repository. Its source was compiled and a
small graph was compared against all eight C++ index files during the preceding
installation audit. Historical whole-genome measurements and annotation-aware
validation are described in the companion repository. No whole-genome build
was repeated for this installation change.

Linux and Windows were not tested in this pass. The documented setup targets
Linux and macOS; the measurements above confirm this macOS ARM64 installation.

## Ubuntu x86-64 validation, 2026-10-02

Native Ubuntu 24.04 x86-64 GitHub-hosted runners tested commit `6df30ea`.
Both Make and CMake built and installed successfully. The Make regression
suite reported **47 passed, 0 failed**; both installations passed the installed
wrappers, Solo matrices, helper lookup and fallback smoke tests.

The first run exposed a missing `<stddef.h>` include in the standalone UMI
header under GCC. Adding that explicit dependency fixed the failure; the
successful rerun is linked below. The workflow now runs on pushes to
`solo/gene-model`, with read-only repository permissions.

[Successful Ubuntu run](https://github.com/iandriver/hisat2/actions/runs/37024820815)

These are correctness and installation checks, not runtime benchmarks. They
exercise HISAT2's small C++ index construction and Solo counting; the separate
Rust builder and whole-human-genome construction were not repeated here.
Windows remains untested in these passes.
