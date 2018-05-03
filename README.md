# MicrostructureNoise

A Julia package for Bayesian volatility estimation in presence of
microstructure noise.

We take a nonparametric Bayesian approach to estimate the volatility in the stochastic differential equation model 

    dXₜ = b(t,Xₜ) dt + s(t) dWₜ, X₀ = x₀, t ∈ [0,T] 

from noisy discrete time observations of its solution

    Yᵢ = X(tᵢ) + Vᵢ,   0 = t₀ < … < tₙ = T, 

and { Vᵢ } denote unobservable stochastic disturbances.

The prior is based on an inverse Gamma Markov chain prior on the squared volatility function and a Kalman filter reconstructing the latent signal.

## Documentation

See [https://mschauer.github.io/MicrostructureNoise.jl/](https://mschauer.github.io/MicrostructureNoise.jl/).

## Installation

Unregistered package. To install, run:

```
Pkg.clone("https://github.com/mschauer/MicrostructureNoise.jl.jl")
```

## Contribute
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/MicrostructureNoise.jl/issues/1) for questions and coordination of the development.

## References

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility learning under microstructure noise. *Preprint*, 2018.
