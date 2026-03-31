# GROMACS fixtures

This directory contains a real and a synthetic/minimal `dhdl.xvg` fixture used to exercise the GROMACS parser paths.

- `minimal_example/dhdl.xvg`: minimal synthetic fixture for focused parser unit and integration coverage
- `lambda_0.xvg`, `lambda_1.xvg`, `lambda_2.xvg`, `lambda_3.xvg`, `lambda_15.xvg`: locally generated GROMACS 2025.2 outputs with a multidimensional `(mass, coul, vdw, bonded, restraint)` schedule, trimmed to 200 data rows each for regression coverage
