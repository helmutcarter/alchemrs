# Getting Started

## Prerequisites

- Rust 1.85 or newer
- AMBER output files 
    - Currently, the parser only support AMBER output files. More are coming. 

The workspace uses Cargo workspaces and is intended to be built and tested with standard Rust tooling.

## Build the workspace

```bash
cargo build --workspace
```

## Run tests

```bash
cargo test --workspace
```

## Build the CLI

```bash
cargo build -p alchemrs-cli --release
```

## Run the top-level examples

```bash
cargo run --example amber_ti -- 300 ./fixtures/amber/acetamide_tiny/\*/acetamide.prod.out

cargo run --example amber_mbar -- 300 ./fixtures/amber/acetamide_tiny/\*/acetamide.prod.out
```

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
