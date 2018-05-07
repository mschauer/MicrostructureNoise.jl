var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MicrostructureNoise.jl-1",
    "page": "Home",
    "title": "MicrostructureNoise.jl",
    "category": "section",
    "text": "MicrostructureNoise is a Julia package for Bayesian volatility estimation in presence of market microstructure noise. The underlying model is the stochastic differential equation $ dX_t=b(t,X_t)\\,dt + s(t)\\,d W_t, \\quad X_0=x_0, \\quad t\\in [0,T] .$The estimation method is minimalistic in its assumptions on the volatility function s, which in particular can be a stochastic process. The process X is latent: observed is its noisy version on a discrete time grid,$ Y_{i}=X_{t_{i}}+V_{i}, \\quad 0<t_1<\\cdots<t_n=T.$Here  V_i  denote unobservable stochastic disturbances, and n is the total number of observations.For data Y_i that are densely spaced in time, the drift function b has little effect on estimation accuracy of the volatility function s, and can be set to zero. This reduces the original model to the linear state space model, and statistical tools developed for the latter can be used to infer the unknown volatility. The posterior inference is performed via the Gibbs sampler, and relies on Kalman filtering ideas to reconstruct unobservable states X(t_i).Essential details of the procedure are as follows: The unknown squared volatility function s^2 is a priori modelled as piecewise constant: Fix an integer mn. Then a unique decomposition n=mN+r with 0leq rm holds, where N=lfloor nmrfloor. Now define bins B_k=t_m(k-1)t_mk), k=1ldotsN-1, and B_N=t_m(N-1)T. The number N of bins is a hyperparameter. Let s be piecewise constant on bins B_k, so that$ s^2=\\sum_{k=1}^{N} \\theta_k \\mathbf{1}_{B_k}.$The coefficients  theta_k  are assigned the inverse Gamma Markov chain prior, which induces smoothing among adjacent pieces of the function s^2. This prior is governed by the smoothing hyperparameter alpha, which in turn is equipped with a hyperprior. The errors V_i are assumed to follow the Gaussian distribution with mean zero and variance eta. The Bayesian model specification is completed by assigning the noise level eta the inverse Gamma prior, and equipping the initial state X_0 with the Gaussian prior. To sample from the joint posterior of the vector theta_k, the noise level eta and the smoothing hyperparameter alpha, the Gibbs sampler is used. In each cycle of the sampler, the unobservable state vector X(t_i) is drawn from its full conditional distribution using the Forward Filtering Backward Simulation algorithm; this employs Kalman filter recursions in the forward pass.Synthetic data examples show that the procedure adapts well to the unknown smoothness of the volatility s.See the referenced article for additional details on prior specification, implementation, and numerical experiments."
},

{
    "location": "index.html#Example-1",
    "page": "Home",
    "title": "Example",
    "category": "section",
    "text": "using MicrostructureNoise, Distributions\n# downloads a large file \nBase.download(\"https://www.truefx.com/dev/data//2015/MARCH-2015/EURUSD-2015-03.zip\",\"./data/EURUSD-2015-03.zip\")\nrun(`unzip ./data/EURUSD-2015-03.zip -d ./data`)\ndat = readcsv(\"./data/EURUSD-2015-03.csv\")\ntimes = map(a -> DateTime(a, \"yyyymmdd H:M:S.s\"), dat[1:10:130260,2])\ntt = Float64[1/1000*(times[i].instant.periods.value-times[1].instant.periods.value) for i in 1:length(times)]\nn = length(tt)-1\nT = tt[end]\ny = Float64.(dat[1:10:130260, 3])\n\nprior = MicrostructureNoise.Prior(\nN = 40,\n\nα1 = 0.0,\nβ1 = 0.0,\n\nαη = 0.3, \nβη = 0.3,\n\nΠα = LogNormal(1., 0.5),\nμ0 = 0.0,\nC0 = 5.0\n)\n\nα = 0.3\nσα = 0.1\ntd, θs, ηs, αs, p = MicrostructureNoise.MCMC(prior, tt, y, α, σα, 1500)\n\nposterior = MicrostructureNoise.posterior_volatility(td, θs)"
},

{
    "location": "index.html#MicrostructureNoise.Prior",
    "page": "Home",
    "title": "MicrostructureNoise.Prior",
    "category": "type",
    "text": "MicrostructureNoise.Prior(; N, α1, β1, αη, βη, Πα, μ0, C0)\nMicrostructureNoise.Prior(; kwargs...)\n\nStruct holding prior distribution parameters. N is the number of bins,  InverseGamma(α1, β1) is the prior of θ[1] on the first bin, the prior on the noise variance η is InverseGamma(αη, βη), the hidden state X_0 at start time is Normal(μ0, C0),  and Πα is a prior Distribution for α,  for example Πα = LogNormal(1., 0.5).\n\nNote: All keyword arguments N, α1, β1, αη, βη, Πα, μ0, C0 are mandatory.\n\nExample:\n\nprior = MicrostructureNoise.Prior(\nN = 40, # number of bins\n\nα1 = 0.0, # prior for the first bin\nβ1 = 0.0,\n\nαη = 0.3, # noise variance prior InverseGamma(αη, βη)\nβη = 0.3,\n\nΠα = LogNormal(1., 0.5),\nμ0 = 0.0,\nC0 = 5.0\n)\n\n\n\n"
},

{
    "location": "index.html#MicrostructureNoise.MCMC",
    "page": "Home",
    "title": "MicrostructureNoise.MCMC",
    "category": "function",
    "text": "MCMC(Π::Union{Prior,Dict}, tt, yy, α0::Float64, σα, iterations; subinds = 1:1:iterations, η0::Float64 = 0.0, printiter = 100) -> td, θ, ηs, αs, pacc\n\nRun the Markov Chain Monte Carlo procedure for iterations iterations, on data (tt, yy), where tt are observation times and yy are observations. α0 is the initial guess for the smoothing parameter α (necessary), η0 is the initial guess for the noise variance (optional), and σα is the stepsize for the random walk proposal for α.\n\nPrints verbose output every printiter iteration.\n\nReturns td, θs, ηs, αs, pacc, td is the time grid of the bin boundaries, ηs, αs are vectors of iterates, possible subsampled at indices subinds, θs is a Matrix with iterates of θ rows. paccα is the acceptance probability for the update step of α.\n\nKeyword args fixalpha, fixeta when set to true allow fixing α and η at their initial values.\n\n\n\n"
},

{
    "location": "index.html#MicrostructureNoise.Posterior",
    "page": "Home",
    "title": "MicrostructureNoise.Posterior",
    "category": "type",
    "text": "struct Posterior\n    post_t # Time grid of the bins\n    post_qlow # Lower boundary of marginal credible band\n    post_median # Posterior median\n    post_qup # Upper boundary of marginal credible band\n    post_mean # Posterior mean of `s^2`\n    post_mean_root # Posterior mean of `s`\n    qu # `qu*100`-% marginal credible band\nend\n\nStruct holding posterior information for squared volatility s^2.\n\n\n\n"
},

{
    "location": "index.html#MicrostructureNoise.posterior_volatility",
    "page": "Home",
    "title": "MicrostructureNoise.posterior_volatility",
    "category": "function",
    "text": "posterior_volatility(td, θs; burnin = size(θs, 2)÷3, qu = 0.90)\n\nComputes the qu*100-% marginal credible band for squared volatility s^2 from θ.\n\nReturns Posterior object with boundaries of the marginal credible band, posterior median and mean of s^2, as well as posterior mean of s.\n\n\n\n"
},

{
    "location": "index.html#MicrostructureNoise.piecewise",
    "page": "Home",
    "title": "MicrostructureNoise.piecewise",
    "category": "function",
    "text": "piecewise(tt, yy, [endtime]) -> tt, xx\n\nIf (tt, yy) is a jump process with piecewise constant paths and jumps  of size yy[i]-y[i-1] at tt[i], piecewise returns coordinates path  for plotting purposes. The second argument allows to choose the right endtime of the last interval.\n\n\n\n"
},

{
    "location": "index.html#Library-1",
    "page": "Home",
    "title": "Library",
    "category": "section",
    "text": "MicrostructureNoise.Prior\nMicrostructureNoise.MCMC\nMicrostructureNoise.Posterior\nMicrostructureNoise.posterior_volatility\nMicrostructureNoise.piecewise"
},

{
    "location": "index.html#Contribute-1",
    "page": "Home",
    "title": "Contribute",
    "category": "section",
    "text": "See issue #1 (Roadmap/Contribution) for questions and coordination of the development."
},

{
    "location": "index.html#References-1",
    "page": "Home",
    "title": "References",
    "category": "section",
    "text": "Shota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility estimation. arxiv:1801.09956, 2018.\nShota Gugushvili, Frank van der Meulen, Moritz Schauer, and Peter Spreij: Nonparametric Bayesian volatility learning under microstructure noise. In preparation."
},

]}
