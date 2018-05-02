function piecewise(tt_, yy, tend = tt[end])
    tt = [tt_[1]]
    n = length(yy)
    append!(tt, repeat(tt_[2:n], inner=2))
    push!(tt, tend)
    tt, repeat(yy, inner=2)
end


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

# Metropolis parameters
function MCMC(Π::Union{Prior,Dict}, tt, y, α0::Float64, σα, iterations; subinds = 1:1:iterations, quc = 0.9, η0::Float64 = 0.0, printiter=100, fixalpha = false, fixeta = false, summaryfile = STDOUT )
    
    N = Π.N

    α1 = Π.α1
    β1 = Π.β1
    αη = Π.αη
    βη = Π.βη
    Πα = Π.Πα 
    μ0 = Π.μ0
    C0 = Π.C0


    Πx0 = Normal(μ0, sqrt(C0))
 

    n = length(y) - 1
    
    N = Π.N
    m = n ÷ N
    


    Πη = InverseGamma(Π.αη, Π.βη)


    # Initialization

    x = copy(y)
    η = η0
    α = α0

    θ = zeros(N)
    ζ = zeros(N - 1)



    Z = zeros(N)
    ii = Vector(N)
    td = zeros(N)
    for k in 1:N
        if k == N
            ii[k] = 1+(k-1)*m:n
        else
            ii[k] = 1+(k-1)*m:(k)*m
        end

        tk = tt[ii[k]]
        td[k] = tk[1]
        Z[k] = sum((x[i+1] - x[i]).^2 ./ (tt[i+1]-tt[i]) for i in ii[k])
        θ[k] = mean(InverseGamma(α1 + length(ii[k])/2, β1 + Z[k]/2))
    end

    acc = 0
    αs = Float64[]
    si = 1
    fixalpha = false
    fixeta = false

    μ = zeros(n+1)
    C = zeros(n+1)
    θs = Any[]
    ηs = Float64[]
    subinds = 1:1:iterations
    samples = zeros(N, length(subinds))

    for iter in 1:iterations
        # update Zk (necessary because x changes)
        for k in 1:N
            Z[k] = sum((x[i+1] - x[i]).^2 ./ (tt[i+1]-tt[i]) for i in ii[k])
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
            push!(αs, α)

        end
        if !fixeta
            # update eta
            assert(length(x) == n + 1)
            z = sum((x[i] - y[i])^2 for i in 2:(n+1))
            η = rand(InverseGamma(αη + n/2, βη + z/2))

            mod(iter, printiter) == 0 && print("\t √η", √(η))


        end
        mod(iter, printiter) == 0 && println()

        # Sample x from Kalman smoothing distribution

        # Forward pass

        C[1] = C0
        μ[1] = μ0
        μi = μ0
        Ci = C0
        for k in 1:N
            iik = ii[k]
            for i in iik # from 1 to n
                wi = θ[k]*(tt[i+1] - tt[i])
                Ki = (Ci + wi)/(Ci + wi + η) # Ci is still C(i-1)
                μi =  μi + Ki*(y[i+1] - μi)

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
                wi1 = θ[k]*(tt[i+1] - tt[i])
                Ci = C[i]
                μi = μ[i]


                Hi = (Ci*wi1)/(Ci + wi1)
                hi = μi + Ci/(Ci + wi1)*(x[i+1] - μi) # x[i+1] was sampled in previous step
                x[i] = rand(Normal(hi, sqrt(Hi)))
            end
        end
        push!(θs, θ)
        push!(ηs, η)
        if iter in subinds
            samples[:, si] = θ
            si += 1
        end
    end
    open(summaryfile, "w") do f
        println(f, "prior alpha", Πα)
        println(f, "n $n nf $nf")
        println(f, "η0 $η0 ")  
        println(f, "N $N m $m")  
        println(f, "α1 $α1 β1 $β1")
        println(f, "αη $αη βη $βη")
        println(f, "x0 ~ N($μ0, $C0)")
        println(f, "σα $σα")
        println(f, "acc α", round(acc/iterations, 3))
        println(f, "fixalpha $fixalpha fixeta $fixeta")
        println(f, "iterations ", iterations)
        println(f, "samples at ", subinds)
        println(f, "$quc % credible bands")
    end

    
    

    samples, ηs, αs
end

struct Posterior
    post_qlow
    post_median
    post_qup
    post_mean
    post_mean_root
end
function compute_posterior_s0(samples; burnin = size(samples, 2)÷3, quc = 0.90)
    p = 1.0 - quc 
    A = view(samples, :, burnin:size(samples, 2))
    post_qup = mapslices(v-> quantile(v, 1 - p/2), A, 2)
    post_mean = mean(A, 2)
    post_mean_root = mean(sqrt.(A), 2)
    post_median = median(A, 2)
    post_qlow = mapslices(v-> quantile(v,  p/2), A, 2)
    Posterior(
        post_qlow,
        post_median,
        post_qup,
        post_mean,
        post_mean_root
    )
end







