#included in bkf.jl

include("delegatemacro.jl")

abstract type AbstractSparseFilter <: KalmanFilter end

DEFAULT_PRODUCT_RADIUS = 1
DEFAULT_RADIUS_FRINGE = 0
DEFAULT_HISTPERLOC = 7


mutable struct FuzzFinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    obs::Union{Observation,Nothing} #the observation which was used to create this. Not used, only bookkeeping.
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
        #[space,phistory]
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
    prevParts = particleMatrix(prev)
    d, n = size(prevParts)
    h = myparams.mh.histPerLoc

    #get new centers
    filt = getbkf(prev)
    means = newCenters(filt, prevParts)


    fuzzes = getNextFuzzes(filt, prevParts, prev.params)
    base = copy(means)

    if myparams.rejuv != 0
      for j = 1:n
          try
            base[:,j] += rand(MvNormal(Matrix(Hermitian(fuzzes[j] * myparams.rejuv))))
                #Matrix(Hermitian( :  ...Work, stupid!
                #/4 : rejuv lightly, but pretend it's full. Not actually correct but meh.
          catch
            base[:,j] += rand(MvNormal(Matrix(Hermitian(fuzzes[j] * myparams.rejuv + Matrix(1e-4I,d,d)))))
          end
      end
    end
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
    lps = Array{Float64,3}(undef,d,n,n)
    for ph = 1:n
        Σ = Matrix(Hermitian(noiseMatrix(filt.f)+fuzzes[ph]))

        for l = 1:d
            hood = prodNeighborhood(l,d,myparams.mh.r)

            try
              localDists[l,ph] = MvNormal(means[hood,ph],Σ[hood,hood])
            catch y
              #
              debug("XXXXX FuzzFinkelParticles error")
              debug(y)
              debug(Σ[hood,hood])
              debug("XXXXXXXX forwardDistribution\n")
              localDists[l,ph] = MvNormal(means[hood,ph],Matrix(1.0I,size(hood)[1],size(hood)[1]))
            end
            #debug("broken",myparams.mh.r,hood,means[hood,ph],localDists[l,ph])
            for pf = 1:n
                lps[l,pf,ph] = logpdf(Normal(means[l,ph],Σ[l,l]),
                                base[l,pf])
            end
        end
    end

    for pf = 1:n
        for l = 1:d
            logForwardDensities[l,pf] = logsumexp(lps[l,pf,:])
            histSampProbs[l,pf] = getSampProbs(lps[l,pf,1:n],pf,myparams.s)
        end
    end
    tip = ParticleSet(filt,
                      n,
                      tipVals,
                      ProbabilityWeights(ones(n)) #dummy value, ignore
                      )
    debugOnce("FinkelParticles", typeof(lps))
    fp = FuzzFinkelParticles(
                tip,
                nothing,
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
    histSampProbs = Array{ProbabilityWeights,2}(undef,d,n)
    localDists = Array{Distribution,2}(undef,0,0)
    logForwardDensities = Array{Float64,2}(undef,0,0)

    lps = Array{Float64,3}(undef,0,0,0)
    tip = ParticleSet(filt,
                      n,
                      tipVals,
                      ProbabilityWeights(ones(n)) #dummy value, ignore
                      )
    debugOnce("FinkelParticles", typeof(lps))
    fp = FuzzFinkelParticles(
                tip,
                nothing,
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

function neighborhoodCenter(fp)
  div(fp.params.mh.r,2) + 1
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
        l::Int64,
        lstem::Nothing)
    pdf(fp.localDists[l,h],fp.tip.particles[neighborhood,i])
end

function probSum(fp::FuzzFinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Vector,
        l::Int64, #optional: location to replace
        lstem::Int64)
    lp = 0.
    state = fp.tip.particles[neighborhood,i]
    state[neighborhoodCenter(fp)] = lstem
    pdf(fp.localDists[l,h],state)
end
