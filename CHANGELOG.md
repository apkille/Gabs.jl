# News

## v1.1.1 - dev

- Add a benchmark suite as a part of the Github workflows.

## v1.1.0 - 2024-11-18

- Add `prob` function for `Generaldyne` type.
- Add features for random generation of Gaussian states, unitaries, and channels.
- **(breaking)** Switched basis for matrices generated by `symplecticform` from block diagonal to canonical.

## v1.0.2 - 2024-11-03

- Remove StaticArrays as a dependency and add as an extension.
- Add `Generaldyne` type with corresponding `output` function.
- Add `randstate` and `randchannel` to calculate random instances of Gabs types.

## v1.0.1 - 2024-10-16

- Add Makie attributes for `heatmap`.

## v1.0.0 - 2024-09-24

- First release.
- Implement `GaussianState`, `GaussianUnitary`, and `GaussianChannel` types.
- Define operations on custom types from `QuantumInterface`, including `ptrace`, `⊗`, and `apply!`. 
- Define predefined functions that create custom types.
- Create Makie extension.
