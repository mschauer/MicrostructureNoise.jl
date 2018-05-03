# MicrostructureNoise

## Overview

MicrostructureNoise is a Julia package for Bayesian volatility estimation in presence of market microstructure noise.

## Installation

Unregistered package. To install, run:

```
Pkg.clone("https://github.com/mschauer/MicrostructureNoise.jl.jl")
```

## Description

MicrostructureNoise estimates the volatility function s of the stochastic differential equation

    dXₜ = b(t,Xₜ) dt + s(t) dWₜ, X₀ = x₀, t ∈ [0,T] 

from noisy observations of its solution

    Yᵢ = X(tᵢ) + Vᵢ,   0 = t₀ < … < tₙ = T, 

where { Vᵢ } denote unobservable stochastic disturbances. The method is minimalistic in its assumptions on the volatility function, which in particular can itself be a stochastic process.

The estimation methodology is intuitive to understand, given that its ingredients are well-known statistical techniques. The posterior inference is performed via the Gibbs sampler, with the Forward Filtering Backward Simulation algorithm used to reconstruct unobservable states X(tᵢ). This relies on the Kalman filter. The unknown squared volatility function is a priori modelled as piecewise constant and is assigned the inverse Gamma Markov chain prior, which induces smoothing among adjacent pieces of the function. The picture below gives an idea of the results obtainable with the method. Depicted is a reconstruction of the volatility function from the synthetic data generated according to the classical Heston stochastic volatility model. Note that next to the point estimate (posterior mean plotted in black), the method conducts automatic uncertainty quantification via the marginal Bayesian credible band (plotted in blue).

<img src="./heston.png" width=600>

## Documentation

See [https://mschauer.github.io/MicrostructureNoise.jl/latest](https://mschauer.github.io/MicrostructureNoise.jl/latest).

## Contribute
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/MicrostructureNoise.jl/issues/1) for questions and coordination of the development.

## References

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility estimation. [arxiv:1801.09956](https://arxiv.org/abs/1801.09956), 2018.

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility learning under microstructure noise. In preparation.
