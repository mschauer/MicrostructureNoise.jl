using MicrostructureNoise, Distributions

using Random, LinearAlgebra

Random.seed!(10)
η = 0.01

n = 1000
tt = 2*[0.0; sort(rand(n-1)); 1.0]
T = tt[end]
σ(t) = 1.0 + sin(t)
# ad hoc euler
x = [0.0; cumsum(σ.(tt[2:end]).*randn(n).*sqrt.(diff(tt)))]
y = x + sqrt(η)*randn(n+1)
prior = MicrostructureNoise.Prior(
N = 20,

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

@test_throws DimensionMismatch MicrostructureNoise.MCMC(prior, tt, y[2:end], α, σα, 10) 


@testset "skipfirst" begin
    @test_logs (:info, "skip observation y[1] at t[1] (skipfirst == true)") MicrostructureNoise.MCMC(prior, tt, y, α, σα, 10; skipfirst = true)
    @test_throws DimensionMismatch MicrostructureNoise.MCMC(prior, tt, y[3:end], α, σα, 10; skipfirst = true) 

    Random.seed!(10)
    p1 = MicrostructureNoise.MCMC(prior, tt, [y[2];y[2:end]], α, σα, 10; skipfirst = true)
    Random.seed!(10)
    p2 = MicrostructureNoise.MCMC(prior, tt, y[2:end], α, σα, 10; skipfirst = true)
    @test p1[1] == p2[1]
    @test p1[2][end] == p2[2][end]
end

posterior = MicrostructureNoise.posterior_volatility(td, θs)

# plot(MicrostructureNoise.piecewise(posterior.post_t[1:end], posterior.post_mean_root[1:end])...)
# plot!(tt,1+sin.(tt), ylim=(0,2))

# Test posterior mean versus truth

@test abs(mean(ηs)-η) < 0.2*η
@test norm(diff(td).*(mean(θs, dims=2) - σ.(td[1:end-1]).^2)) < 1