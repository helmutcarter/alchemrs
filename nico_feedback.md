Knowing nothing about the domain, this looks well-organized and quite clear. I would suggest updating the README examples with some specific example commands on particular data files so that the prospective user can get a sample of how it's intended to work. I tried to do so with random prod.out files, and got errors, such as this one:

```
alchemrs on  main is 📦 v0.1.0 via 🦀 v1.94.0 took 2s
❯ cargo run -p alchemrs-cli --release -- ti \
  --temperature 300 \
  --method trapezoidal \
  --remove-burnin 125 \
  fixtures/amber/acetamide_tiny/0.7/acetamide.prod.out

    Finished `release` profile [optimized] target(s) in 0.02s
     Running `target/release/alchemrs-cli ti --temperature 300 --method trapezoidal --remove-burnin 125 fixtures/amber/acetamide_tiny/0.7/acetamide.prod.out`
Error: "remove-burnin exceeds series length"
```
 or this one:

```
❯ cargo run -p alchemrs-cli --release -- ti \
  --temperature 300 \
  --method trapezoidal \
  --remove-burnin 0 \
  fixtures/amber/acetamide_tiny/1.0/acetamide.prod.out

    Finished `release` profile [optimized] target(s) in 0.02s
     Running `target/release/alchemrs-cli ti --temperature 300 --method trapezoidal --remove-burnin 0 fixtures/amber/acetamide_tiny/1.0/acetamide.prod.out`
Error: InvalidShape { expected: 2, found: 1 }
```

I expect the reasons why these fail would be more obvious to the target audience, but it's rather unclear to me, even looking through the source code. It's not clear if I'm passing in the wrong type of file, if the specific combination of cli options are incompatible, or what.

I had assumed that the .prod.out files were the right ones to pass in as arguments to the cli based on the README, though looking through them, this may not be the case. I also tried the `ti_0.0_1.0.delta_f_sigma.txt` file, which (I think?) fixed that error but resulted in a different error further down the pipeline.

```
❯ cargo run -p alchemrs-cli --release -- bar \
  --temperature 300 \
  --method false-position \
  --decorrelate \
  fixtures/amber/acetamide_tiny/ti_0.0_1.0.delta_f_sigma.txt

Error: Parse("failed to locate temp0 in AMBER output")
```

Looking at some of the possible error points, it would seem that the InvalidShape error 

```rust
if values.len() < 2 {
    return Err(CoreError::InvalidShape {
        expected: 2,
        found: values.len(),
    });
}
```

is triggered when a 1-d array is of length < 2, ie either 0 or 1, but would not trigger if the length was more than 2. Coming from a Python background, my instinct when seeing an error like `InvalidShape { expected: 2, found: 1 }` is that it's a dimension mismatch error, and that it's expecting 2-dimensional data, but is getting 1-dimensional data. In this case, the data is meant to be 1-dimensional, and the error is not in the 'shape' as commonly construed, but in the length, and 2 is not the only valid lenght, but it must be at least 2. An error like "Provided data must have at least two data points" could help be more informative, assuming that is, in fact, the intended return.

The construction sites of the InvalidShape error could do to be more informative, especially since it's not at all clear to the user what the cause of the error was, at what point in the pipeline, or what they can do to change this.

Also, the *prod.out files have 'file generated for tiny fixtures' repeated multiple times at the top, which I assume is not the intended output. 

## Linter warnings and errors

There are also some spots in the code where the linter is catching some small complaints, which can mostly be fixed automatically, with a couple exceptions. Consider the code from extract_u_nk in alchemrs-parse:

```rust
for sample_idx in 0..n_samples {
    let mut ok = true;
    for state_idx in 0..n_states {
        if !mbar_energies[state_idx][sample_idx].is_finite() {
            ok = false;
            break;
        }
        if ok {
            valid_indices.push(sample_idx);
        }
    }
}
```

Can be rewritten as

```rust
(0..n_samples).for_each(|sample_idx| {
    if (0..n_states).all(|state_idx| {
        mbar_energies[state_idx][sample_idx].is_finite()
    }) {
        valid_indices.push(sample_idx);
    }
});
```

for concision and clarity, which then also fixes some linter errors in 1.94.0.

Asides the linter warnings, there are also some errors which do not prevent compilation, but do still look quite odd, such as these two lines in `alchemrs-cli/main.rs`:

```rust
let delta = result.values()[0 * n + (n - 1)];

let delta = result.values()[(n - 1) * n + 0];
```

The linter complains that, in the former case, the result of the `0 * n` expression will always evaluate to 0, and in the latter case, that the `+ 0` addition has no effect. I assume you left these here as a reminder for something else later, but it escapes me what that could be. Also highlighting these two examples because they both look quite similar, and clippy volunteers to 'fix' both automatically, but the corrected expressions are different, `n-1` and `n(n-1)` respectively. It's unclear which is intended, or if they're intended to be the same or different.

# Conclusion

Overall, I think this could do with another manual pass to look for more obtuse code like the double loop above, and to clarify the error output to be more communicative to the end user.

