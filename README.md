[![Build Status](https://travis-ci.org/mschauer/MicrostructureNoise.jl.svg?branch=master)](https://travis-ci.org/mschauer/MicrostructureNoise.jl)
[![Coverage Status](https://coveralls.io/repos/github/mschauer/MicrostructureNoise.jl/badge.svg?branch=master)](https://coveralls.io/github/mschauer/MicrostructureNoise.jl?branch=master)
[![codecov.io](http://codecov.io/github/mschauer/MicrostructureNoise.jl/coverage.svg?branch=master)](http://codecov.io/github/mschauer/MicrostructureNoise.jl?branch=master)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://mschauer.github.io/MicrostructureNoise.jl/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1241011.svg)](https://doi.org/10.5281/zenodo.1241011)


# MicrostructureNoise

## Overview

MicrostructureNoise is a [Julia](https://github.com/JuliaLang/julia) package for Bayesian volatility estimation in presence of market microstructure noise.

## Installation

To install, run:

```
Pkg.add("MicrostructureNoise")
```

## Description

MicrostructureNoise estimates the volatility function <a href="https://www.codecogs.com/eqnedit.php?latex=s" target="_blank"><img src="https://latex.codecogs.com/svg.latex?s" title="s" /></a> of the stochastic differential equation

<a href="https://www.codecogs.com/eqnedit.php?latex=dX_t&space;=&space;b(t,X_t)&space;dt&space;&plus;&space;s(t)&space;dW_t,&space;\quad&space;X_0&space;=&space;x_0,&space;\quad&space;t&space;\in&space;[0,T]" target="_blank"><img src="https://latex.codecogs.com/svg.latex?dX_t&space;=&space;b(t,X_t)&space;dt&space;&plus;&space;s(t)&space;dW_t,&space;\quad&space;X_0&space;=&space;x_0,&space;\quad&space;t&space;\in&space;[0,T]" title="dX_t = b(t,X_t) dt + s(t) dW_t, \quad X_0 = x_0, \quad t \in [0,T]" /></a>

from noisy observations of its solution

<a href="https://www.codecogs.com/eqnedit.php?latex=Y_i&space;=&space;X(t_i)&space;&plus;&space;V_i,&space;\quad&space;0&space;<&space;t_1&space;<&space;\ldots&space;<&space;t_n&space;=&space;T," target="_blank"><img src="https://latex.codecogs.com/gif.latex?Y_i&space;=&space;X(t_i)&space;&plus;&space;V_i,&space;\quad&space;0&space;<&space;t_1&space;<&space;\ldots&space;<&space;t_n&space;=&space;T," title="Y_i = X(t_i) + V_i, \quad 0 < t_1 < \ldots < t_n = T," /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\{V_i\}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\{V_i\}" title="\{V_i\}" /></a> denote unobservable stochastic disturbances. The method is minimalistic in its assumptions on the volatility function, which in particular can itself be a stochastic process.

The estimation methodology is intuitive to understand, given that its ingredients are well-known statistical techniques. The posterior inference is performed via the Gibbs sampler, with the Forward Filtering Backward Simulation algorithm used to reconstruct unobservable states <a href="https://www.codecogs.com/eqnedit.php?latex=X(t_i)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?X(t_i)" title="X(t_i)" /></a>. This relies on the Kalman filter. The unknown squared volatility function is a priori modelled as piecewise constant and is assigned the inverse Gamma Markov chain prior, which induces smoothing among adjacent pieces of the function. The picture below gives an idea of the results obtainable with the method. Depicted is a reconstruction of the volatility function from the synthetic data generated according to the classical Heston stochastic volatility model (the unobserved true volatility curve is plotted in red). Note that next to the point estimate (posterior mean plotted in black), the method conducts automatic uncertainty quantification via the marginal Bayesian credible band (plotted in blue).

<img src="./heston.png" width=600>

When <a href="https://www.codecogs.com/eqnedit.php?latex=X(t_i)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?X(t_i)" title="X(t_i)" /></a> is observed without noise, an option (`fixeta`) allows to perfom inference as described in the reference "Nonparametric Bayesian volatility estimation".

## Documentation

See [https://mschauer.github.io/MicrostructureNoise.jl/latest](https://mschauer.github.io/MicrostructureNoise.jl/latest).

## Contribute
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/MicrostructureNoise.jl/issues/1) for questions and coordination of the development.

## References

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility estimation. [arxiv:1801.09956](https://arxiv.org/abs/1801.09956), 2018.

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility learning under microstructure noise. [arxiv:1805.05606](https://arxiv.org/abs/1805.05606), 2018.

* A. T. Cemgil and O. Dikmen: Conjugate gamma Markov random fields for modelling nonstationary sources. In ICA 2007, 7th International Conference on Independent Component Analysis and Signal Separation, September 2007.
