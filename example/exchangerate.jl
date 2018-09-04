
using MicrostructureNoise, Distributions
using Plots

using Statistics, DelimitedFiles, Dates
# uncomment if you do not mind to create this large file 
# Base.download("https://www.truefx.com/dev/data//2015/MARCH-2015/EURUSD-2015-03.zip","./data/EURUSD-2015-03.zip")
# run(`unzip ./data/EURUSD-2015-03.zip -d ./data`)
dat = readdlm("./data/EURUSD-2015-03.csv", ',')
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

plot(MicrostructureNoise.piecewise(posterior.post_t, posterior.post_mean_root)..., label="Volatility estimate")
