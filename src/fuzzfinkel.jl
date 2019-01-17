#included in bkf.jl

include("delegatemacro.jl")

abstract type AbstractSparseFilter <: KalmanFilter end

DEFAULT_PRODUCT_RADIUS = 1
DEFAULT_RADIUS_FRINGE = 0
DEFAULT_HISTPERLOC = 7


mutable struct FuzzFinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    base::Array{T,2} #This is the raw 1-step progression from last time.
                #As with all similar arrays, size is [d,n]; that is, first index is space and second is particle.
        #[space, particle]
    prev::AbstractFinkel
    historyTerms::Array{Int64,3} #tells which histories each particle was checked against
        #[space, particle, sample]
    stem::Array{Int64,2} #tells which base each tip comes from
        #[space, particle]
    ws::Vector{ProbabilityWeights} #selection probabilities at each point; p(y_l|z^i_l)
        #[space][particle]
    logForwardDensities::Array{Float64,2} #forward densities at each point; p(z_l|x^{1..M}_l)
        #sum of lps over phistory
        #[space, particle]
    lps::Array{Float64,3} #log probabilities from hist to fut; size [d, n, n]
        #[space, pfuture, phistory]
        #assuming all other loci are at mean
    histSampProbs::Array{ProbabilityWeights,2}
        #[space, pfuture][phistory]
    means::Array{T,2} #base, progressed, without adding noise.
        #[space, particle]
    prevProbs::Array{Union{Float64,Nothing},2} #sum over neighborhood history of local probability - avoid double calculation
    localDists::Array{Distribution,2} #for calculating probs
        #[space, phistory]
    totalProb::Array{Float64,2} #convergence diagnostic
    params::P
    numMhAccepts::Int64


end


function FuzzFinkelParticles(prev::AbstractFinkel)
    myparams = fparams(prev)
    FuzzFinkelParticles(prev,myparams)
end

function FuzzFinkelParticles(prev::AbstractFinkel, nosteps::Nothing)
    myparams = fparams(prev)
    FuzzFinkelParticles(prev,myparams,nosteps)
end

function FuzzFinkelParticles(prev::AbstractFinkel,
                         myparams::FinkelParams)
    d, n = size(particleMatrix(prev))
    h = myparams.mh.histPerLoc

    #get new centers
    filt = getbkf(prev)
    prevParts = particleMatrix(prev)
    means = newCenters(filt, prevParts)


    base = ap(prev.tip).particles #TODO: should use means instead of recalculating them
    fuzzes = getNextFuzzes(filt, prevParts, prev.params)
    tipVals = copy(base)


    historyTerms = zeros(Int64,d,n,h)
    stem = zeros(Int64,d,n)
    for j = 1:n
        for l = 1:d
            stem[l,j] = j
            for η = 1:h
                historyTerms[l,j,η] = j
            end
        end
    end

    histSampProbs = Array{ProbabilityWeights,2}(undef,d,n)
    localDists = Array{Distribution,2}(undef,d,n)
    logForwardDensities = Array{Float64,2}(undef,d,n)
    for ph = 1:n
        for l = 1:d #Ideally we could somehow call forwardDistribution just once but meh
            localDists[l,ph] = forwardDistribution(filt.f,
                                    means[:,ph],
                                    fuzzes[ph],
                                    l)
        end
    end

    lps = Array{Float64,3}(undef,d,n,n)
    for pf = 1:n
        for l = 1:d

            for ph = 1:n
                lps[l,pf,ph] = logpdf(localDists[l,ph],
                                base[l:l,pf])
            end
            logForwardDensities[l,pf] = logsumexp(lps[l,pf,:])
            histSampProbs[l,pf] = getSampProbs(lps[l,pf,1:n],pf,myparams.s)
        end
    end
    tip = ParticleSet(filt,
                      n,
                      tipVals,
                      ProbabilityWeights(ones(n)) #dummy value, ignore
                      )
    if params.rejuv
      tip=rejuvenate(tip)
    end
    debugOnce("FinkelParticles", typeof(lps))
    fp = FuzzFinkelParticles(
                tip,
                base, #base
                prev,
                historyTerms, #historyTerms
                stem,
                Vector{ProbabilityWeights}(), #empty weights
                logForwardDensities,
                lps,
                histSampProbs,
                means, #means
                Array{Union{Nothing, Float64},2}(nothing,d,n), #prevProbs
                localDists, #localDists
                zeros(d,n), #totalProb
                myparams,
                0 #numMhAccepts
                )
    #getDists!(fp, d, n)
    #calcPrevProbs!(fp, d, n)
    fp
end

function getNextFuzzes(filt, prevParts, params)
  d, n = size(prevParts)
  Σ = cov(prevParts;dims=2)
  Σinv = inv(Σ)
  μ = mean(prevParts;dims=2)
  localVariance = Array{Float64,2}(undef,n,d) #particle, locus
  basePrecision = Vector{Float64}(undef,d)
  diffs = prevParts .- μ
  for l in 1:d
    hood = prodNeighborhood(l,d,params.mh.r)
    hooddet = det(Σ[hood,hood])
    dethat = prod(Σ[i,i] for i in hood)
    basePrecision[l] = Σ[l,l] * (hooddet/dethat)^(1/params.mh.r)
    for p in 1:n
      dist = MvNormal(Σ[hood,hood])
      logRelDens = logpdf(dist,diffs[hood,p]) - logpdf(dist,zeros(params.mh.r))
      localVariance[p,l] = 1/(n-params.mh.r)/(1+exp(logRelDens))
    end
  end
  [propagateUncertainty(filt.f,prevParts[:,p:p],
      inv(Σinv + Diagonal([basePrecision[l] / params.overlap *
        prod(exp(localVariance[p,l]/params.mh.r^2)
        for λ in prodNeighborhood(l,d,params.mh.r))
      for l in 1:d])))
    for p in 1:n]
end

function FuzzFinkelParticles(prev::AbstractFinkel,
                         myparams::FinkelParams,
                         nosteps::Nothing)
    d, n = size(particleMatrix(prev))
    h = myparams.mh.histPerLoc
    base = ap(prev.tip).particles
    tipVals = copy(base)


    historyTerms = zeros(Int64,0,0,0)
    stem = zeros(Int64,d,n)

    filt = getbkf(prev)
    means = newCenters(filt, particleMatrix(prev))
    histSampProbs = Array{ProbabilityWeights,2}(d,n)
    localDists = Array{Distribution,2}(0,0)
    logForwardDensities = Array{Float64,2}(0,0)

    lps = Array{Float64,3}(0,0,0)
    tip = ParticleSet(filt,
                      n,
                      tipVals,
                      ProbabilityWeights(ones(n)) #dummy value, ignore
                      )
    debugOnce("FinkelParticles", typeof(lps))
    fp = FinkelParticles(
                tip,
                base, #base
                prev,
                historyTerms, #historyTerms
                stem,
                Vector{ProbabilityWeights}(), #empty weights
                logForwardDensities,
                lps,
                histSampProbs,
                means, #means
                fill(Nullable{Float64}(),d,n), #prevProbs
                localDists, #localDists
                zeros(d,n), #totalProb
                myparams,
                0 #numMhAccepts
                )
    fp
end

"""
    probSum(fp,i,h,
            neighborhood,

            l=-1,lstem=0)

    multiply likelihood from h to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(fp::FuzzFinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Vector,
        l::Nothing,
        lstem::Nothing)
    pdf(fp.localDists[h],fp.tip.particles[neighborhood,i])
end

function probSum(fp::FuzzFinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Vector,
        l::Int64, #optional: location to replace
        lstem::Int64)
    lp = 0.
    state = fp.tip.particles[neighborhood,i]
    state[findfirst(x -> x==l,neighborhood)] = lstem
    pdf(fp.localDists[h],state)
end
