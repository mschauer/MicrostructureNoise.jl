
"""
    MicrostructureNoise.Prior(; N, α1, β1, αη, βη, Πα, μ0, C0)
    MicrostructureNoise.Prior(; kwargs...)

Struct holding prior distribution parameters.
`N` is the number of bins, 
`InverseGamma(α1, β1)` is the prior of `θ[1]` on the first bin,
the prior on the noise variance `η` is `InverseGamma(αη, βη)`,
the hidden state ``X_0`` at start time is `Normal(μ0, C0)`, 
and `Πα` is a prior `Distribution` for `α`, 
for example `Πα = LogNormal(1., 0.5)`.

Note: All keyword arguments `N, α1, β1, αη, βη, Πα, μ0, C0`
are mandatory.


Example:
```
prior = MicrostructureNoise.Prior(
N = 40, # number of bins

α1 = 0.0, # prior for the first bin
β1 = 0.0,

αη = 0.3, # noise variance prior InverseGamma(αη, βη)
βη = 0.3,

Πα = LogNormal(1., 0.5),
μ0 = 0.0,
C0 = 5.0
)
```
"""
struct Prior
    N

    α1
    β1

    αη
    βη

    Πα 

    μ0
    C0

    function Prior(;Π_...)
        Π = Dict(Π_)
        return new(
            Π[:N],
            Π[:α1],
            Π[:β1],
            Π[:αη],
            Π[:βη],
            Π[:Πα],
            Π[:μ0],
            Π[:C0]
            )
    end
end

"""
    MCMC(Π::Union{Prior,Dict}, t, y, α0::Float64, σα, iterations; 
        subinds = 1:1:iterations, η0::Float64 = 0.0, printiter = 100,
        fixalpha = false, fixeta = false, skipfirst = false) -> td, θ, ηs, αs, pacc

Run the Markov Chain Monte Carlo procedure for `iterations` iterations,
on data `(t, y)`, where `t` are observation times and `y` are observations.
`α0` is the initial guess for the smoothing parameter `α` (necessary),
`η0` is the initial guess for the noise variance (optional),
and `σα` is the stepsize for the random walk proposal for `α`.

Prints verbose output every `printiter` iteration.

Returns `td, θs, ηs, αs, pacc`,
`td` is the time grid of the bin boundaries,
`ηs`, `αs` are vectors of iterates,
possible subsampled at indices `subinds`,
`θs` is a Matrix with iterates of `θ` rows.
`paccα` is the acceptance probability for the update step of `α`.

`y[i]` is the observation at `t[i]`.

If `skipfirst = true` and `t` and `y` are of equal length,
the observation `y[1]` (corresponding to `t[1]`) is ignored.

If `skipfirst = true` and `length(t) = length(y) + 1`, 
`y[i]` is the observation at `t[i + 1]`.

Keyword args `fixalpha`, `fixeta` when set to `true` allow fixing
`α` and `η` at their initial values. 
"""
function MCMC(Π::Union{Prior,Dict}, t, y, α0::Float64, σα, iterations; subinds = 1:1:iterations, η0::Float64 = 0.0, printiter = 100, 
    fixalpha = false, fixeta = false, skipfirst = false)
    
    N = Π.N

    α1 = Π.α1
    β1 = Π.β1
    αη = Π.αη
    βη = Π.βη
    Πα = Π.Πα 
    μ0 = Π.μ0
    C0 = Π.C0


    Πx0 = Normal(μ0, sqrt(C0))
 
    if skipfirst
        if length(t) == length(y)
            shift = 0
            @info "skip observation y[1] at t[1] (skipfirst == true)"
        elseif length(t) == length(y) + 1
            shift = 1
        else
            throw(DimensionMismatch("expected length(t) or length(t) - 1 observations (skipfirst == true)"))   
        end      
    else
        length(t) != length(y) && throw(DimensionMismatch("expected length(t) observations"))   
        shift = 0
    end

    n = length(t) - 1 # number of increments  
    
    N = Π.N
    m = n ÷ N
    


    Πη = InverseGamma(Π.αη, Π.βη)


    # Initialization
    if shift == 1
        x = [y[1]; y]
    else
        x = copy(y)
    end
    η = η0
    α = α0

    θ = zeros(N)
    ζ = zeros(N - 1)



    Z = zeros(N)
    ii = Vector(undef, N) # vector of start indices of increments
    td = zeros(N+1)
    for k in 1:N
        if k == N
            ii[k] = 1+(k-1)*m:n # sic!
        else
            ii[k] = 1+(k-1)*m:(k)*m
        end

        tk = t[ii[k]]
        td[k] = tk[1]
        Z[k] = sum((x[i+1] - x[i]).^2 ./ (t[i+1]-t[i]) for i in ii[k])
        θ[k] = mean(InverseGamma(α1 + length(ii[k])/2, β1 + Z[k]/2))
    end
    td[end] = t[end]

    acc = 0
    αs = Float64[]
    si = 1
 
    μ = zeros(n+1)
    C = zeros(n+1)
    θs = Any[]
    ηs = Float64[]

    samples = zeros(N, length(subinds))

    for iter in 1:iterations
        # update Zk (necessary because x changes)
        if !(η == 0.0 && fixeta)
            for k in 1:N
                Z[k] = sum((x[i+1] - x[i]).^2 ./ (t[i+1]-t[i]) for i in ii[k])
            end
        end

        # sample chain
        for k in 1:N-1
            ζ[k] = rand(InverseGamma(α + α, (α/θ[k] + α/θ[k+1])))
        end
        for k in 2:N-1
            θ[k] = rand(InverseGamma(α + α + m/2, (α/ζ[k-1] + α/ζ[k] + Z[k]/2)))
        end
        θ[1] = rand(InverseGamma(α1 + α + m/2, β1 + α/ζ[1] + Z[1]/2))
        θ[N] = rand(InverseGamma(α + length(ii[N])/2, α/ζ[N-1] + Z[N]/2))
   
        if !fixalpha
            # update parameter alpha using Wilkinson II
            α˚ = α + σα*randn()
            while α˚ < eps()
                α˚ = α + σα*randn()
            end

            lq = logpdf(Πα, α)
            lq += (2*(N-1))*(α*log(α) - lgamma(α))
            s = sum(log(θ[k-1]*θ[k]*ζ[k-1]*ζ[k-1]) + (1/θ[k-1] + 1/θ[k])/ζ[k-1]  for k = 2:N)
            lq += -α*s

            lq˚ = logpdf(Πα, α˚)
            lq˚ += (2*(N-1))*(α˚*log(α˚) - lgamma(α˚))
            lq˚ += -α˚*s

            mod(iter, printiter) == 0 && print("$iter \t α ", α˚)
            if rand() < exp(lq˚ - lq)*cdf(Normal(0, σα), α)/cdf(Normal(0, σα), α˚) # correct for support
                acc = acc + 1
                α = α˚
                mod(iter, printiter) == 0 && print("✓")
            end
            

        end
        if !fixeta
            # update eta
            z = sum((x[i] - y[i - shift])^2 for i in (1 + skipfirst):(n+1))
            η = rand(InverseGamma(αη + n/2, βη + z/2))

            mod(iter, printiter) == 0 && print("\t √η", √(η))


        end
        mod(iter, printiter) == 0 && println()

        # Sample x from Kalman smoothing distribution

        # Forward pass
        if η == 0.0 && fixeta
            # do nothing
        else 
            if skipfirst # ignore observation y[1]
                C[1] = C0
                μ[1] = μ0
                μi = μ0
                Ci = C0
            else
                wi = 0.0 
                Ki = C0/(C0 + η)
                μi =  μ0 + Ki*(y[1] - μ0) # shift = 0
   
                Ci = Ki*η
                C[1] = Ci
                μ[1] = μi
            end
            for k in 1:N
                iik = ii[k]
                for i in iik # from 1 to n
                    wi = θ[k]*(t[i+1] - t[i])
                    Ki = (Ci + wi)/(Ci + wi + η) # Ci is still C(i-1)
                    μi =  μi + Ki*(y[i + 1 - shift] - μi)

                    Ci = Ki*η
                    C[i+1] = Ci # C0 is at C[1], Cn is at C[n+1] etc.
                    μ[i+1] = μi
                end
            end
            hi = μi
            Hi = Ci

            # Backward pass
            x[end] = rand(Normal(hi, sqrt(Hi)))

            for k in N:-1:1
                iik = ii[k]
                for i in iik[end]:-1:iik[1] # n to 1
                    wi1 = θ[k]*(t[i+1] - t[i])
                    Ci = C[i]
                    μi = μ[i]


                    Hi = (Ci*wi1)/(Ci + wi1)
                    hi = μi + Ci/(Ci + wi1)*(x[i+1] - μi) # x[i+1] was sampled in previous step
                    x[i] = rand(Normal(hi, sqrt(Hi)))
                end
            end
        end
        if iter in subinds
            push!(αs, α)
            push!(ηs, η)
            samples[:, si] = θ
            si += 1
        end
    end

    td, samples, ηs, αs, round(acc/iterations, digits=3)
end

"""
```
struct Posterior
    post_t # Time grid of the bins
    post_qlow # Lower boundary of marginal credible band
    post_median # Posterior median
    post_qup # Upper boundary of marginal credible band
    post_mean # Posterior mean of `s^2`
    post_mean_root # Posterior mean of `s`
    qu # `qu*100`-% marginal credible band
end
```

Struct holding posterior information for squared volatility `s^2`.
"""
struct Posterior
    post_t
    post_qlow
    post_median
    post_qup
    post_mean
    post_mean_root
    qu
end

"""
    posterior_volatility(td, θs; burnin = size(θs, 2)÷3, qu = 0.90)

Computes the `qu*100`-% marginal credible band for squared volatility `s^2` from `θ`.

Returns `Posterior` object with boundaries of the marginal credible band,
posterior median and mean of `s^2`, as well as posterior mean of `s`.
"""
function posterior_volatility(td, samples; burnin = size(samples, 2)÷3, qu = 0.90)
    p = 1.0 - qu 
    A = view(samples, :, burnin:size(samples, 2))
    post_qup = mapslices(v-> quantile(v, 1 - p/2), A, dims=2)
    post_mean = mean(A, dims=2)
    post_mean_root = mean(sqrt.(A), dims=2)
    post_median = median(A, dims=2)
    post_qlow = mapslices(v-> quantile(v,  p/2), A, dims=2)
    Posterior(
        td,
        post_qlow[:],
        post_median[:],
        post_qup[:],
        post_mean[:],
        post_mean_root[:],
        qu
    )
end

"""
    piecewise(t, y, [endtime]) -> t, xx

If `(t, y)` is a jump process with piecewise constant paths and jumps 
of size `y[i]-y[i-1]` at `t[i]`, piecewise returns coordinates path 
for plotting purposes. The second argument
allows to choose the right endtime of the last interval.
"""
function piecewise(t_, y, tend = t_[end])
    t = [t_[1]]
    n = length(y)
    append!(t, repeat(t_[2:n], inner=2))
    push!(t, tend)
    t, repeat(y, inner=2)
end





