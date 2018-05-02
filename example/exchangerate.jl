
using MicrostructureNoise, Distributions
if !isdefined(:dat)
    dat = readcsv("./data/EURUSD-2015-03.csv")
end
times = map(a -> DateTime(a, "yyyymmdd H:M:S.s"), dat[1:10:130260,2])
tt = Float64[1/1000*(times[i].instant.periods.value-times[1].instant.periods.value) for i in 1:length(times)]
n = length(tt)-1
T = tt[end]
y = Float64.(dat[1:10:130260, 3])

Π = MicrostructureNoise.Prior(
N = 40,
α = 0.1,
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
MicrostructureNoise.MCMC(Π, tt, y, α, σα, 10000)

