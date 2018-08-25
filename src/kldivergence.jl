function kl2(dist::MvNormal,
            samps::Array{Float64,2},window = 3)
    p = params(dist)
    μ1 = p[1]
    μ2 = vec(mean(samps,2))
    d = length(μ1)
    Σ1 = p[2] * eye(d)
    Σ2 = cov(samps,2)
    kl2(μ1,Σ1,μ2,Σ2, window)
end

function kl2(μ1,Σ1,μ2,Σ2, window = 3)

    d = length(μ1)
    nwind = div(d,window)
    subdivs = zeros(nwind)
    subtr = zeros(nwind)
    subdif = zeros(nwind)
    sublog = zeros(nwind)
    for w = 1:nwind
        r = ((w-1) * window + 1):(w * window)
        diff = (μ2[r] - μ1[r])
        subtr[w] = (trace(inv(Σ2[r,r]) * Σ1[r,r]) - window) / 2
        subdif[w] = diff' * Σ2[r,r] * diff / 2
        d2, d1 = det(Σ2[r,r]), det(Σ1[r,r])
        if (d2 <= 0)
            print("kl2 non-positive determinant", d1, " ", d2, " ")
            print(Σ1[r,r])
            print(Σ2[r,r])
            d2 = d1
        end
        sublog[w] = (log(d2/d1)) / 2
        subdivs[w] = subtr[w] + subdif[w] + sublog[w]
    end
    #print("KL parts:",mean(subtr),",",mean(subdif),",",mean(sublog),",",mean(subdivs))
    [mean(subdivs),mean(subtr),mean(subdif),mean(sublog)]
end


function sqerr(truth::Array{Float64,1},
            samps::Array{Float64,2})
    nsamps = size(samps)[2]
    mean(
        mean((truth - samps[:,i]) .^ 2) for i = 1:nsamps
        )
end

function sqerr(truth::Array{Float64,1},
    f::Union{ParticleStep,FrankenStep})
    sqerr(truth,f.p)
end

function sqerr(truth::Array{Float64,1},
    f::ParticleSet)
    samps = f.particles

    nsamps = size(samps)[2]
    mean(
        [mean((truth - samps[:,i]) .^ 2) for i = 1:nsamps], f.weights
    )
end

function sqerr(truth::Array{Float64,1},
    f::FrankenSet)

    w = nlocs(f)
    mu = zeros(w)
    sig = zeros((w,w))

    lim = size(f.weights,1)
    submeans=zeros(lim)
    for n in 1:lim
        toN = n*f.hoodSize
        fromN = toN - f.hoodSize + 1
        # print(n," ",size(f.weights[n])," ",f.n," ")
        # print(mean((truth[fromN:toN,i] - f.particles[fromN:toN,1]) .^ 2)," ")
        # print()
        submeans[n] = mean(
            [mean((truth[fromN:toN] - f.particles[fromN:toN,i]) .^ 2) for i = 1:f.n], f.weights[n]
        )
    end
    mean(submeans)
end

#http://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html
function logsumexp_batch(X)
    alpha = maximum(X)  # Find maximum value in X
    log(sum(exp(X-alpha))) + alpha
end

#http://www.nowozin.net/sebastian/blog/streaming-log-sum-exp-computation.html
function logsumexp_stream(X)
    alpha = -Inf
    r = 0.0
    for x = X
        if x <= alpha
            r += exp(x - alpha)
        else
            r *= exp(alpha - x)
            r += 1.0
            alpha = x
        end
    end
    log(r) + alpha
end

logsumexp = logsumexp_batch

function meanvarlocs(f::ParticleStep, r::Range)
    locs = vec(sum(f.p.particles[r,:],1))
    (mean(locs, f.p.weights),
        var(locs, f.p.weights))
end

function meanvarlocs(f::FinkelParticles, r::Range)
    meanvarlocs(f.tip,r)
end


function meanvarlocs(f::Union{BasicKalmanFilter, FrankenStep}, r::Range)
    p=musig(f)
    meanvarlocs(p...,r)
end

function meanvarlocs(f::Observation, r::Range)
    (sum(f.y[r]),0.)
end

function meanvarlocs(f::ParticleSet, r::Range)
    if f.n == 1
        (sum(f.particles[r,1]),0.)
    else
        locs = vec(sum(f.particles[r,:],1))
        (mean(locs),
            var(locs))
    end
end

function meanvarlocs(f::FrankenSet, r::Range)
    (sum(f.particles[r,1]),0.)
end

function meanvarlocs(μ,Σ,r::Range)
    d = length(μ)
    v = zeros(d)
    v[r] = 1.

    (sum(v .* μ), v⋅(Σ*v))
end


function musig(f::FrankenSet, lim=999)

    w = nlocs(f)
    mu = zeros(w)
    sig = zeros((w,w))

    for i in 1:min(size(f.weights,1),lim)
        toI = i*f.hoodSize
        fromI = toI - f.hoodSize + 1

        mu[fromI:toI] = mean(f.particles[fromI:toI,:],f.weights[i],2)
        sig[fromI:toI,fromI:toI] = cov(f.particles[fromI:toI,:],f.weights[i],2, corrected=false)
    end
    (mu,sig)
end

function musig(f::FrankenStep, lim=999)
    musig(f.p, lim)
end

function musig(f::BasicKalmanFilter)
    musig(toDistribution(f))
end

function musig(f::MvNormal)
    p = params(f)
    (p[1],full(p[2]))
end

function musig(f::ParticleStep)
    (mean(f.p.particles,f.p.weights),
      cov(f.p.particles,f.p.weights,corrected=false))
end
