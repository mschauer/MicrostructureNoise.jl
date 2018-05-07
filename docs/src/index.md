# MicrostructureNoise.jl

`MicrostructureNoise` is a Julia package for Bayesian volatility estimation in presence of market microstructure noise. The underlying model is the stochastic differential equation 

$ dX_t=b(t,X_t)\,dt + s(t)\,d W_t, \quad X_0=x_0, \quad t\in [0,T] .$

The estimation method is minimalistic in its assumptions on the volatility function $s$, which in particular can be a stochastic process. The process $X$ is latent: observed is its noisy version on a discrete time grid,

$ Y_{i}=X_{t_{i}}+V_{i}, \quad 0<t_1<\cdots<t_n=T.$

Here $\{ V_i \}$ denote unobservable stochastic disturbances, and $n$ is the total number of observations.

For data $\{Y_i\}$ that are densely spaced in time, the drift function $b$ has little effect on estimation accuracy of the volatility function $s$, and can be set to zero. This reduces the original model to the linear state space model, and statistical tools developed for the latter can be used to infer the unknown volatility. The posterior inference is performed via the Gibbs sampler, and relies on Kalman filtering ideas to reconstruct unobservable states $\{X(t_i)\}$.

Essential details of the procedure are as follows: The unknown squared volatility function $s^2$ is a priori modelled as piecewise constant: Fix an integer $m<n$. Then a unique decomposition $n=mN+r$ with $0\leq r<m$ holds, where $N=\lfloor {n}/{m}\rfloor$. Now define bins $B_k=[t_{m(k-1)},t_{mk})$, $k=1,\ldots,N-1$, and $B_N=[t_{m(N-1)},T]$.
The number $N$ of bins is a hyperparameter. Let $s$ be piecewise constant on bins $B_k$, so that

$ s^2=\sum_{k=1}^{N} \theta_k \mathbf{1}_{B_k}.$

The coefficients $\{ \theta_k \}$ are assigned the inverse Gamma Markov chain prior, which induces smoothing among adjacent pieces of the function $s^2$. This prior is governed by the smoothing hyperparameter $\alpha$, which in turn is equipped with a hyperprior. The errors $\{V_i\}$ are assumed to follow the Gaussian distribution with mean zero and variance $\eta$. The Bayesian model specification is completed by assigning the noise level $\eta$ the inverse Gamma prior, and equipping the initial state $X_0$ with the Gaussian prior. To sample from the joint posterior of the vector $\{\theta_k\}$, the noise level $\eta$ and the smoothing hyperparameter $\alpha$, the Gibbs sampler is used. In each cycle of the sampler, the unobservable state vector $\{X(t_i)\}$ is drawn from its full conditional distribution using the Forward Filtering Backward Simulation algorithm; this employs Kalman filter recursions in the forward pass.

Synthetic data examples show that the procedure adapts well to the unknown smoothness of the volatility $s$.

See the referenced article for additional details on prior specification, implementation, and numerical experiments.

## Example

```
using MicrostructureNoise, Distributions
# downloads a large file 
Base.download("https://www.truefx.com/dev/data//2015/MARCH-2015/EURUSD-2015-03.zip","./data/EURUSD-2015-03.zip")
run(`unzip ./data/EURUSD-2015-03.zip -d ./data`)
dat = readcsv("./data/EURUSD-2015-03.csv")
times = map(a -> DateTime(a, "yyyymmdd H:M:S.s"), dat[1:10:130260,2])
tt = Float64[1/1000*(times[i].instant.periods.value-times[1].instant.periods.value) for i in 1:length(times)]
n = length(tt)-1
T = tt[end]
y = Float64.(dat[1:10:130260, 3])

prior = MicrostructureNoise.Prior(
N = 40,

α1 = 0.0,
β1 = 0.0,

αη = 0.3, 
βη = 0.3,

Πα = LogNormal(1., 0.5),
μ0 = 0.0,
C0 = 5.0
)

α = 0.3
σα = 0.1
td, θs, ηs, αs, p = MicrostructureNoise.MCMC(prior, tt, y, α, σα, 1500)

posterior = MicrostructureNoise.posterior_volatility(td, θs)
```

## Library

```@docs
MicrostructureNoise.Prior
MicrostructureNoise.MCMC
MicrostructureNoise.Posterior
MicrostructureNoise.posterior_volatility
MicrostructureNoise.piecewise
```

## Contribute
See [issue #1 (Roadmap/Contribution)](https://github.com/mschauer/MicrostructureNoise.jl/issues/1) for questions and coordination of the development.

## References

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility estimation. [arxiv:1801.09956](https://arxiv.org/abs/1801.09956), 2018.

* Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility learning under microstructure noise. In preparation.
