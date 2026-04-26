# Getting Started

The fastest way to use `alchemrs` is through the CLI. The Rust API is documented separately for embedding and custom workflows.

## Prerequisites

- Rust 1.85 or newer
- AMBER output files or GROMACS `dhdl.xvg` files
- Python 3.9 or newer plus `maturin` if you want the local Python package

The repository is a normal Cargo workspace and can be built and tested with standard Rust tooling.

## Build the package

```bash
cargo build
```

## Build the CLI binary

```bash
cargo build --release
```

Then run:

```bash
./target/release/alchemrs --help
```

## Try a CLI workflow

```bash
./target/release/alchemrs ti \
  --temperature 300 \
  ./fixtures/amber/acetamide_tiny/*/acetamide.prod.out
```

## Run tests

```bash
cargo test
```

## Run the Rust examples

```bash
cargo run --example amber_ti -- 300 ./fixtures/amber/acetamide_tiny/*/acetamide.prod.out

cargo run --example amber_mbar -- 300 ./fixtures/amber/acetamide_tiny/*/acetamide.prod.out

cargo run --example openmm_u_kln_mbar
```

## Run the Python package locally

Build and install the Python extension into the active environment:

```bash
maturin develop
```

Then run the Python tests and example scripts from the repository root:

```bash
python -m pytest python/tests -q
python python/examples/amber_fixture_analysis.py
python python/examples/atm_analysis.py
```

The OpenMM toy-system examples are also available when `openmm` is installed:

```bash
python python/examples/openmm_u_kln_mbar.py
python python/examples/openmm_nes.py
python python/examples/openmm_atm.py
```

See [Python and OpenMM](python.md) for the Python package layout, helper APIs,
and the current local-usage workflow.

## Build this book

This documentation uses `mdBook`.

Install it once:

```bash
cargo install mdbook
```

Preview the docs locally:

```bash
mdbook serve docs
```

Build static HTML:

```bash
mdbook build docs
```
