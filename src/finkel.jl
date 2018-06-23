#included in bkf.jl

include("delegatemacro.jl")

abstract type AbstractSparseFilter <: KalmanFilter end

DEFAULT_PRODUCT_RADIUS = 2
DEFAULT_HISTPERLOC = 2

type SparseKF <: AbstractSparseFilter
    f::BasicKalmanFilter
    neighbors::Vector{Vector{Int16}}
end






@delegate SparseKF.f [ covs, toDistribution, noiseDistribution, obsNoiseDistribution ]

function SparseKF(f::BasicKalmanFilter) #assumes neighborliness is symmetric
    neighbors = Vector{Vector{Int16}}()
    m = f.f.a
    for j = 1:size(m,2)
        nset = Vector{Int16}()
        for i = 1:size(m,1)
            if m[i,j] != 0
                push!(nset, i)
            end
        end
        push!(neighbors, nset)
    end
    SparseKF(f,neighbors)
end

abstract type AbstractFinkel <: AbstractParticleFilter end

type FinkelToe{T,F<:KalmanFilter} <: AbstractFinkel
    tip::ParticleSet{T,F}
end

function particleMatrix(fp::AbstractFinkel)
    fp.tip.particles
end

function getbkf(fp::AbstractFinkel)
    fp.tip.filter
end

abstract type SampleType end #how to sample particle histories

type SampleUniform <: SampleType
end

type SampleLog <: SampleType
    inflectionPoint::Float64 #Distance below max log prob where probability slope bends
    factor::Float64 #how much probability slope bends
end

abstract type MhType end #how to calculate MH acceptance probability

type MhLocal <: MhType #Restricted product, unrestricted sum
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
end

type MhSampled <: MhType #Neighborhood product, sum using locally sampled history
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples
end

type MhMultilocus <: MhType #Variable neighborhood product, sum using histories sampled for each center
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples per locus
end

type MhMultisampled <: MhType #Fixed neighborhood product, sum using histories sampled for each locus
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples per locus
end

type MhCompromise <: MhType #Fixed neighborhood product, sum using histories sampled for each locus in core subset
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples per locus
    rSub::Int64 #neighborhood size for history samples
end


type FinkelParams{S<:SampleType,MH<:MhType}
    s::S
    mh::MH
end

function fparams(basedOn::Any)
    FinkelParams(SampleUniform(),MhSampled(DEFAULT_PRODUCT_RADIUS,DEFAULT_HISTPERLOC))
end

type FinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    base::Array{T,2} #This is the raw 1-step progression from last time.
                #As with all similar arrays, size is [d,n]; that is, first index is space and second is particle.
        #[space, particle]
    prev::AbstractFinkel
    historyTerms::Array{Int64,3} #tells which histories each particle was checked against
        #[space, sample, particle]
    stem::Array{Int64,2} #tells which base each tip comes from
        #[space, particle]
    ws::Vector{ProbabilityWeights} #selection probabilities at each point; p(y_l|x^i_l)
        #[space][particle]
    lps::Array{Float64,3} #log probabilities from hist to fut; size [d, n, n]
        #[space, pfuture, phistory]
    histSampProbs::Array{ProbabilityWeights,2}
        #[space, pfuture][phistory]
    means::Array{T,2} #base, progressed, without adding noise.
        #[space, particle]
    prevProbs::Array{Nullable{Float64},2} #sum over neighborhood history of local probability - avoid double calculation
    localDists::Array{Distribution,2} #for calculating probs
        #[space, phistory]
    totalProb::Array{Float64,2} #convergence diagnostic
    params::P
    numMhAccepts::Int64
end

function fparams(fp::FinkelParticles)
    fp.params
end


function FinkelParticles(prev::AbstractFinkel,
                         )
    d, n = size(particleMatrix(prev))
    myparams = fparams(prev)
    h = myparams.mh.histPerLoc
    print(d," dn ",n,"\n")
    base = ap(prev.tip).particles
    tipVals = copy(base)


    historyTerms = zeros(Int64,d,h,n)
    stem = zeros(Int64,d,n)
    for j = 1:n
        for l = 1:d
            stem[l,j] = j
            for η = 1:h
                historyTerms[l,η,j] = j
            end
        end
    end

    filt = getbkf(prev)
    means = filt.f.a * particleMatrix(prev)
    histSampProbs = Array{ProbabilityWeights,2}(d,n)
    localDists = Array{Distribution,2}(d,n)
    for ph = 1:n
        for l = 1:d #Ideally we could somehow call forwardDistribution just once but meh
            localDists[l,ph] = forwardDistribution(getbkf(prev).f,
                                    means[ph,],
                                    l)
        end
    end
    lps = Array{Float64,3}(d,n,n)
    for pf = 1:n
        for l = 1:d
            for ph = 1:n
                lps[l,pf,ph] = logpdf(localDists[l,ph],
                                base[l,pf])
            end
            histSampProbs[l,pf] = getSampProbs(lps[l,pf,1:n],myparams.s)
        end
    end
    tip = ParticleSet(filt,
                      n,
                      tipVals,
                      ProbabilityWeights(ones(n)) #dummy value, ignore
                      )
    fp = FinkelParticles(
                tip,
                base, #base
                prev,
                historyTerms, #historyTerms
                stem,
                Vector{ProbabilityWeights}(), #empty weights
                lps,
                histSampProbs,
                means, #means
                Array{Nullable{Float64},2}(d,n), #prevProbs
                localDists, #localDists
                zeros(d,n), #totalProb
                myparams,
                0 #numMhAccepts
                )
    getDists!(fp, d, n)
    calcPrevProbs!(fp, d, n)
    fp
end

function getSampProbs(logpdfs,
                    s::SampleUniform)
    ProbabilityWeights(ones(logpdfs))
end

function particles(fp::FinkelParticles)
    fp.tip.particles
end

function reweight!(fp::FinkelParticles, y::Observation)
    fp.ws = Vector{ProbabilityWeights}()
    diffs = fp.base - fp.tip.filter.z.h * repeat(y.y,outer=[1,fp.tip.n])  # fp.tip.f.z.h should probably be eye ?
    vars = diag(fp.tip.filter.z.r) #assumes fp.tip.f.z.r is diagonal
    for i = 1:size(fp.base,1)
        wvec = Vector{Float64}(fp.tip.n)
        for j = 1:fp.tip.n
            wvec[j] = exp(-diffs[i,j]^2 / vars[i] / 2)
        end
        push!(fp.ws,ProbabilityWeights(wvec))
    end
end

function replant!(fp::FinkelParticles)
    d = size(fp.base,1)
    for i in 1:fp.tip.n
      for l in 1:d
          p = sample(fp.ws[l])

          fp.stem[l,i] = p
          for h in 1:fparams(fp).mh.histPerLoc
              fp.historyTerms[l,h,i] = p
          end
          fp.tip.particles[l,i] = fp.base[l,p]
      end
    end
end

function getDists!(fp::FinkelParticles, d, n)
    for i in 1:n
        x = Vector(fp.means[:,i])
        for l in 1:d
            fp.localDists[l,i] = forwardDistribution(getbkf(fp).f,x,max(l-1,1):min(l+1,d))
        end
    end
end

function histTerms(l,
                    p,
                    fp)
    sample(1:fp.tip.n,
            fp.histSampProbs[l,p],
            fp.params.mh.histPerLoc)
end

"""
    probSum(fp,i,h,
            neighborhood,

            l=-1,lstem=0)

    multiply likelihood from h to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(fp::FinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Range,
        l::Void,
        lstem::Void)
    lp = 0.
    for λ in neighborhood
        lp += fp.lps[λ,fp.stem[λ,i],h]
    end
    exp(-lp)
end

function probSum(fp::FinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Range,
        l::Int64, #optional: location to replace
        lstem::Int64)
    lp = 0.
    for λ in neighborhood
        if λ != l  ## TODO: optimize; check at start, not each step
            lp += fp.lps[λ,fp.stem[λ,i],h]
        else
            lp += fp.lps[λ,lstem,h]
        end
    end
    exp(-lp)
end

"""
    probSum...

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(fp::FinkelParticles{},
    i::Int64,#current particle
    d::Int64,
    ln::Int64, #location for neighborhood center
    lh::Int64, #location for history samples
    lr::T = nothing, #location to replace
    lstem::T = nothing,
    lhist::Void = nothing) where T <: Union{Void,Int64}

    myneighborhood = prodNeighborhood( ln, fp, d, fp.params.mh.r)
    prob = 0.
    for h in fp.historyTerms[lh,:,i]
        prob += probSum(fp, i, h, myneighborhood,
            lr, lstem
            )
    end
    prob
end

function probSum(fp::FinkelParticles{},
    i::Int64,#current particle
    d::Int64,
    ln::Int64, #location for neighborhood center
    lh::Int64, #location for history samples
    lr::T, #location to replace
    lstem::T,
    lhist::Vector{Int64}) where T <: Union{Void,Int64}

    myneighborhood = prodNeighborhood( ln, fp, d, fp.params.mh.r)
    prob = 0.
    if lr == lh
        lhist = Vector(fp.historyTerms[lh,:,i])
    end
    for h in lhist
        prob += probSum(fp, i, h, myneighborhood,
            lr, lstem
            )
    end
    prob
end

"""
    probSum for MhSampled

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::FinkelParticles{T,F,FinkelParams{S,MhSampled}},
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Void,Int64},
                                XT2 <: Union{Void,Vector{Int64}},
                                T, F, S}

    if typeof(lstem) == Void
        probSum(fp,i,d,l,l,lstem,lstem,lhist)
    else
        probSum(fp,i,d,l,l,l,lstem,lhist)
    end
end


function clearOldProbs!(
    fp::FinkelParticles{T,F,FinkelParams{S,MhSampled}},
    i::Int64,#current particle
    l::Int64,
    d::Int64) where {T, F, S}

    #do nothing - mhsampled means it's just local
end


function calcPrevProbs!(fp::FinkelParticles, d, n)

    for i in 1:n
        for l in 1:d
            fp.prevProbs[l,i] = Nullable{Float64}()#Nullable(probSum(fp,l,i,d))
        end
    end
end

function hNeighborsOf(l::Int64, #neighborhood over which we use historyTerms
                    fp::FinkelParticles,
                    d::Int64 #dimension/number of loci
                    )
    if l == 1
        Set(2)
    elseif l == d
        Set(l-1)
    else
        Set(l-1,l+1)
    end
end

function prodNeighborhood(l::Int64, #neighborhood over which we use historyTerms
                    fp::FinkelParticles,
                    d::Int64, #dimension/number of loci
                    w=2)
    max(l-w,1):min(l+w,d)
end


function mcmc!(fp::FinkelParticles,i::Int64,steps::Int64)
    d = size(fp.base,1)
    for s in 1:steps
        order = randperm(d)
        for l in order
            p = sample(fp.ws[l])
            if fp.base[l,p] != fp.tip.particles[l,i]
                oldProbNull = fp.prevProbs[l,i]
                if isnull(oldProbNull)
                    oldProb = probSum(i, d, l, fp)
                    fp.prevProbs[l,i] = Nullable(oldProb)
                else
                    oldProb = get(oldProbNull)
                end
                newHistoryTerms = histTerms(l, p, fp)
                newProb = probSum(i, d, l, fp, p, newHistoryTerms)
                if (newProb < oldProb) && (newProb < rand() * oldProb)
                    fp.totalProb[l,i] += newProb / oldProb
                    continue #M-H rejection
                end
                #M-H accepted

                fp.numMhAccepts += 1
                fp.totalProb[l,i] += 1
                fp.stem[l,i] = p
                fp.historyTerms[l,:,i] = newHistoryTerms

                fp.tip.particles[l,i] = fp.base[l,p]
                fp.prevProbs[l,i] = newProb
                clearOldProbs!(fp,i,l,d)
            end
        end
    end
end

function FinkelParticles(prev::AbstractFinkel, y::Observation, nIter=15, debug=true)
    fp = FinkelParticles(prev)
    reweight!(fp, y)
    replant!(fp)
    for i in 1:fp.tip.n
        mcmc!(fp,i,nIter)
        if debug
            print("Ran particle ", i, "; mean tp = ", mean(fp.totalProb[:,1:i]), "\n")
        end
    end
    fp
end
