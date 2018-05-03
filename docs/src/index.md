# MicrostructureNoise.jl

A Julia package for Bayesian volatility estimation in presence of
microstructure noise.

We take a nonparametric Bayesian approach to estimate the volatility in the stochastic differential equation model 

$ dX_t=b(t,X_t)\,dt + s(t)\,d W_t, \quad X_0=x_0, \quad t\in[0,T] $

from noisy discrete time observations of its solution

$ Y_{i}=X_{t_{i}}+V_{i}, \quad 0=t_0<t_1<\cdots<t_n=T, $

and $\{ V_i \}$ denote unobservable stochastic disturbances.


The prior is based on an inverse Gamma Markov chain prior on the squared volatility function and a Kalman filter reconstructing the latent signal.


For the measurement error variance $\eta_v$, we use a (standard) prior specification $\eta_v \sim IG(\alpha_v,\beta_v)$. Fix an integer $m<n$. Then we have a unique decomposition $n=mN+r$  with $0\leq r<m$, where $N=\lfloor {n}/{m}\rfloor$. 

Now define bins $B_k=[t_{m(k-1)},t_{mk})$, $k=1,\ldots,N-1$, and $B_N=[t_{m(N-1)},T]$.

The number $N$ of bins is a hyperparameter of our prior. Let $s$ be piecewise constant on bins $B_k$, so that

$ s^2=\sum_{k=1}^{N} \theta_k \mathbf{1}_{B_k}.$

We use a Markov chain monte carlo procedure to sample from the joint posterior of the vector `\theta`, the noise level `\eta` and the smoothness hyperparameter `\alpha`.

The procedure incorporating a smoothness hyperparameter shows in a simulation study shows good adaptivity to the unkown smoothness of the volatility.

See the referenced article for details of prior specification, implementation and numerical experiments.

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
θs, ηs, αs, p = MicrostructureNoise.MCMC(prior, tt, y, α, σα, 1500)

posterior = MicrostructureNoise.posterior_volatility(θs)
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