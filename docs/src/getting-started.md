# Getting Started

## Prerequisites

- Rust 1.85 or newer
- AMBER output files for the current parsing workflow

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
cargo run --example amber_ti -- 300 path/to/lambda0.out path/to/lambda1.out

cargo run --example amber_mbar -- 300 path/to/lambda0.out path/to/lambda1.out path/to/lambda2.out
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
