# Gabs.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://apkille.github.io/Gabs.jl/stable)
[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://apkille.github.io/Gabs.jl/dev)
[![Build status (Github Actions)](https://github.com/apkille/Gabs.jl/workflows/CI/badge.svg)](https://github.com/apkille/Gabs.jl/actions)
[![codecov](https://codecov.io/github/apkille/Gabs.jl/graph/badge.svg?token=JWMOD4FY6P)](https://codecov.io/github/apkille/Gabs.jl)
[![JET static analysis](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Gabs.jl is a numerical tooling package for simulating **Gaussian quantum information**.

Gaussian states and operators have the convenient property that they can be
characterized by low-dimensional matrices in the phase space representation.
Thus, a large class of continuous variable quantum information can be efficiently
simulated on a classical computer, lending to applications in quantum cryptography, quantum machine learning, integrated quantum photonics, and more. Gabs.jl provides a high-level [Julia](https://julialang.org) interface for straightforward and high performance implementations of Gaussian quantum systems.