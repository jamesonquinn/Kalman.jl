#included in bkf.jl

#include("delegatemacro.jl")
#using Distributed

abstract type AbstractSparseFilter <: KalmanFilter end

DEFAULT_PRODUCT_RADIUS = 4
DEFAULT_RADIUS_FRINGE = 0
DEFAULT_HISTPERLOC = 7

mutable struct SparseKF <: AbstractSparseFilter
    f::BasicKalmanFilter
    neighbors::Vector{Vector{Int16}}
end

#@mydelegate SparseKF.f [ covs, toDistribution, noiseDistribution, obsNoiseDistribution ]
function covs(skf::SparseKF,x...)
    covs(skf.f,x...)
end



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

abstract type SampleType end #how to sample particle histories

mutable struct SampleUniform <: SampleType
end

mutable struct SampleLog <: SampleType
    inflectionPoint::Float64 #Distance below max log prob where probability slope bends
    factor::Float64 #how much probability slope bends
end

function SampleLog()
    SampleLog(5., 5.) #defaluts
end

mutable struct SampleLogM <: SampleType
    inflectionPoint::Float64 #Distance below max log prob where probability slope bends
    factor::Float64 #how much probability slope bends
    minSpan::Float64 #minimum distance from max to min
end

function SampleLogM(inflectionPoint::Float64,factor::Float64)
    SampleLog(inflectionPoint,factor,inflectionPoint*1.5)
end

abstract type MhType end #how to calculate MH acceptance probability

mutable struct MhLocal <: MhType #Restricted product, unrestricted sum
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
end

mutable struct MhSampled <: MhType #Neighborhood product, sum using locally sampled history
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples
end

function MhSampled(histPerLoc)
    MhSampled(DEFAULT_PRODUCT_RADIUS, histPerLoc)
end

mutable struct MhMultilocus <: MhType #Variable neighborhood product, sum using histories sampled for each center
    neighbors::Int64 #neighborhood size for product.
    histPerLoc::Int64 #history samples per locus
end

mutable struct MhMultisampled <: MhType #Fixed neighborhood product, sum using histories sampled for each locus
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples per locus
end

mutable struct MhCompromise <: MhType #Fixed neighborhood product, sum using histories sampled for each locus in core subset
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples per locus
    rSub::Int64 #neighborhood size for history samples
end

function MhCompromise(histPerLoc) #actually, total histories per attempt, rounded to nearest
    MhCompromise(DEFAULT_PRODUCT_RADIUS, round(Int64,histPerLoc/
                        (DEFAULT_PRODUCT_RADIUS - DEFAULT_RADIUS_FRINGE)),
                DEFAULT_PRODUCT_RADIUS - DEFAULT_RADIUS_FRINGE)
end

mutable struct FinkelParams{S<:SampleType,MH<:MhType}
    s::S
    mh::MH
    useForward::Float64
    overlap::Float64 #overlap factor for fuzz
    algo::Type
    rejuv::Float64
end


abstract type AbstractFinkel <: AbstractParticleFilter end

mutable struct FinkelToe{T,F<:KalmanFilter} <: AbstractFinkel
    tip::ParticleSet{T,F}
    params::FinkelParams
end
#
# function params(ft::FinkelToe)
#   ft.params
# end

function FinkelToe(model::KalmanFilter, n, params)
    FinkelToe(ParticleSet(model,n),params)
end

function particleMatrix(fp::AbstractFinkel)
    fp.tip.particles
end

function getbkf(fp::AbstractFinkel)
    fp.tip.filter
end

#uniform, sampled
function fparams(histPerLoc::Int64 = DEFAULT_HISTPERLOC,
            radius::Int64 = DEFAULT_PRODUCT_RADIUS,
            useForward::Float64=1.,
            overlap::Float64=2.,
            algo::Type=FinkelParticles,
            rejuv=.25)
    FinkelParams(SampleUniform(),MhSampled(radius, histPerLoc),
                useForward, overlap, algo,rejuv)
end

#uniform, compromise
function fparams(
            radius::Int64 ,
            histPerLoc::Int64 ,
            rSub::Int64,
            useForward::Float64=1.,
            overlap::Float64=2.,
            algo::Type=FinkelParticles,
            rejuv=.25
            )
    FinkelParams(SampleUniform(),
                MhCompromise(radius, histPerLoc, rSub),
                useForward, overlap, algo, rejuv)
end

#log, sampled
function fparams(histPerLoc::Int64,
            inflectionPoint::Float64,
            factor::Float64,
            useForward::Float64=1.,
            overlap::Float64=2.,
            algo::Type=FinkelParticles,
            rejuv=.25)
    FinkelParams(SampleLog(inflectionPoint,factor),MhSampled(DEFAULT_PRODUCT_RADIUS, histPerLoc),
                useForward, overlap, algo,
                rejuv)
end

#log, compromise
function fparams(inflectionPoint::Float64,
            factor::Float64,

            radius::Int64 ,
            histPerLoc::Int64 ,
            rSub::Int64,
            useForward::Float64=1.,
            overlap::Float64=2.,
            algo::Type=FinkelParticles,
            rejuv=.25
            )
    FinkelParams(SampleLog(inflectionPoint,factor),
                MhCompromise(radius, histPerLoc, rSub),
                useForward, overlap, algo,
                rejuv)
end



#uniform, sampled is default
function fparams(basedOn::Any)
    fparams()
end

mutable struct FinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    obs::Union{Observation,Nothing} #the observation which was used to create this. Not used, only bookkeeping.
    base::Array{T,2} #This is the raw 1-step progression from last time.
                #As with all similar arrays, size is [d,n]; that is, first index is space and second is particle.
        #[space, particle]
    prev::AbstractFinkel
    historyTerms::SharedArray{Int64,3} #tells which histories each particle was checked against
        #[space, particle, sample]
    stem::SharedArray{Int64,2} #tells which base each tip comes from
        #[space, particle]
    ws::Vector{ProbabilityWeights} #selection probabilities at each point; p(y_l|z^i_l)
        #[space][particle]
    logForwardDensities::Array{Float64,2} #forward densities at each point; p(z_l|x^{1..M}_l)
        #sum of lps
        #[space, particle]
    lps::Array{Float64,3} #log probabilities from hist to fut; size [d, n, n]
        #[space, pfuture, phistory]
    histSampProbs::Array{ProbabilityWeights,2}
        #[space, pfuture][phistory]
    means::Array{T,2} #base, progressed, without adding noise.
        #[space, particle]
    prevProbs::SharedArray{Union{Float64,Nothing},2} #sum over neighborhood history of local probability - avoid double calculation
    #localDists::Array{Distribution,2} #for calculating probs
        #[space, phistory]
    totalProb::SharedArray{Float64,2} #convergence diagnostic
    params::P
    numMhAccepts::Int64
end

function fparams(fp::AbstractFinkel)
    fp.params
end

function FinkelParticles(prev::AbstractFinkel)
    myparams = fparams(prev)
    FinkelParticles(prev,myparams)
end

function FinkelParticles(prev::AbstractFinkel, nosteps::Nothing)
    myparams = fparams(prev)
    FinkelParticles(prev,myparams,nosteps)
end

function FinkelParticles(prev::AbstractFinkel,
                         myparams::FinkelParams)
    d, n = size(particleMatrix(prev))
    h = myparams.mh.histPerLoc

    #get new centers
    filt = getbkf(prev)
    prevParts = particleMatrix(prev)
    means = newCenters(filt, prevParts)


    base = ap(prev.tip).particles #TODO: should use means instead of recalculating them
    fuzz = getNextFuzz(filt, prevParts, prev.params.mh.r)
    tipVals = copy(base)


    historyTerms = SharedArray(zeros(Int64,d,n,h))
    stem = SharedArray(zeros(Int64,d,n))
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
                                    fuzz,
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
    debugOnce("FinkelParticles", typeof(lps))
    fp = FinkelParticles(
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
                SharedArray{Union{Nothing, Float64},2}(nothing,d,n), #prevProbs
                #localDists, #localDists
                SharedArray(zeros(d,n)), #totalProb
                myparams,
                0 #numMhAccepts
                )
    #getDists!(fp, d, n)
    #calcPrevProbs!(fp, d, n)
    fp
end


function FinkelParticles(prev::AbstractFinkel,
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
    #localDists = Array{Distribution,2}(0,0)
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
                #localDists, #localDists
                zeros(d,n), #totalProb
                myparams,
                0 #numMhAccepts
                )
    fp
end

function getSampProbs(logpdfs,
                    pf::Int64,
                    s::SampleUniform)
    ProbabilityWeights(fill(1.,size(logpdfs)))
end

function getSampProbs(lpdfs,
                    pf::Int64,
                    s::SampleLog)
    logpdfs = copy(lpdfs)
    mn, mx = extrema(logpdfs)
    mn -= 1
    #logpdfs[pf] = mn #avoid (over-juicing low-density futures by choosing the past they're conditional on)
    if (mx - mn) < s.inflectionPoint
        result = logpdfs .- mn
    else
        v = logpdfs .- mn
        inflec = mx - mn - s.inflectionPoint
        for i = 1:length(v)
            if v[i] > inflec
                v[i] += (v[i] - inflec) * s.factor
            end
        end
        result = v
    end
    replace_nan!(result)
    ProbabilityWeights(result)
end

function particles(fp::AbstractFinkel)
    fp.tip.particles
end

function reweight!(fp::AbstractFinkel, y::Observation)
    d = size(fp.base,1)
    fp.ws = Vector{ProbabilityWeights}(undef, d)
    diffs = fp.base - fp.tip.filter.z.h * repeat(y.y,outer=[1,fp.tip.n])  # fp.tip.f.z.h should probably be eye ?
    vars = diag(fp.tip.filter.z.r) #assumes fp.tip.f.z.h is eye and ...r is diagonal
    for l = 1:d
        if false
            forwardProb = zeros(d)
            for ph = 1:fp.tip.n
                forwardProb += exp.(fp.lps[l,:,ph])
            end
        else
            forwardProb = 1. #ones(d)
        end
        wvec = exp.(-diffs[l,:].^2 ./ vars[l] / 2 ) ./ forwardProb # ./ if forwardProb is vector

        replace_nan!(wvec)
        fp.ws[l] = ProbabilityWeights(wvec)
    end
    fp.obs=y
end

function replant!(fp::AbstractFinkel)
    d = size(fp.base,1)
    for i in 1:fp.tip.n
      for l in 1:d
          p = sample(fp.ws[l])

          fp.stem[l,i] = p

          if size(fp.historyTerms)[1] > 0
              fp.historyTerms[l,i,:] = histTerms(l, p, fp)
          end
          fp.tip.particles[l,i] = fp.base[l,p]
      end
    end
end
#
# function getDists!(fp::FinkelParticles, d, n)
#     for i in 1:n
#         x = Vector(fp.means[:,i])
#         for l in 1:d
#             fp.localDists[l,i] = forwardDistribution(getbkf(fp).f,x,max(l-1,1):min(l+1,d))
#         end
#     end
# end

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
        neighborhood::Vector,
        l::Nothing,
        lstem::Nothing)
    lp = 0.
    for λ in neighborhood
        lp += fp.lps[λ,fp.stem[λ,i],h]
    end
    exp(lp)
end

function probSum(fp::FinkelParticles,
        i::Int64, #current particle
        h::Int64, #history
        neighborhood::Vector,
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
    exp(lp)
end

function getSampProb(fp::AbstractFinkel,
        i::Int64,#current particle
        l::Int64, #location for neighborhood center
        lstem::Int64)

    fp.histSampProbs[l,i][lstem] / sum(fp.histSampProbs[l,i])
end

function getSampProb(fp::AbstractFinkel,
        i::Int64,#current particle
        ln::Int64, #location for neighborhood center
        lr::Nothing = nothing, #location to replace
        lstem::Nothing = nothing)
    getSampProb(fp,i,ln,fp.stem[ln,i])
end


function getSampProb(fp::FinkelParticles{T,F,FinkelParams{S,MhSampled}},
        i::Int64,#current particle
        ln::Int64, #location for neighborhood center
        myneighborhood, #location for history samples
        lr::Int64, #location to replace
        lstem::Int64) where {T,F,S}

    if lr==ln
        getSampProb(fp,i,ln,lstem)
    else
        getSampProb(fp,i,ln)
    end
end

function getSampProb(fp::AbstractFinkel     ,#{T,F,FinkelParams{S,MhMultisampled}},
        i::Int64,#current particle
        ln::Int64, #location for neighborhood center
        myneighborhood::Vector,## was UnitRange{Int64},
        lr::Union{Int64,Nothing}, #location to replace
        lstem::Union{Int64,Nothing}) #where {T,F,S}

    result = 0.
    for l in myneighborhood
        if lr==l
            result += getSampProb(fp,i,ln,lstem)
        else
            result += getSampProb(fp,i,ln)
        end
    end
    result
end

"""
    probSum...

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(fp::AbstractFinkel,
    i::Int64,#current particle
    d::Int64,
    locn::Int64, #location for neighborhood center
    lh::Int64, #location for history samples
    lr::Int64, #location to replace
    lstem::T = nothing,
    lhist::Nothing = nothing) where T <: Union{Nothing,Int64}

    myneighborhood = prodNeighborhood( locn, d, fp.params.mh.r)
    prob = 0.
    for h in fp.historyTerms[lh,i,:]
        prob += probSum(fp, i, h, myneighborhood,
            lr, lstem
            ) / getSampProb(fp, i, locn, myneighborhood, lr, lstem)
    end
    prob
end

function probSum(fp::AbstractFinkel,
    i::Int64,#current particle
    d::Int64,
    ln::Int64, #location for neighborhood center
    lh::Int64, #location for history samples
    lr::T, #location to replace
    lstem::T,
    lhist::Vector{Int64}) where T <: Union{Nothing,Int64}

    myneighborhood = prodNeighborhood( ln, d, fp.params.mh.r)
    prob = 0.
    if lr == lh
        lhist = Vector(fp.historyTerms[lh,i,:])
    end
    for h in lhist
        prob += probSum(fp, i, h, myneighborhood,
            lr, lstem
            )
    end
    prob
end

function probSum(
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::AbstractFinkel,
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Nothing,Int64}, #lhist is the newly proposed history samples
                                XT2 <: Union{Nothing,Vector{Int64}},
                                T, F, S}
    probSum(fp.params,i,d,l,fp,lstem,lhist)
end

"""
    probSum for MhSampled

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(
    p::FinkelParams{S,MhSampled},
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::AbstractFinkel,
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Nothing,Int64}, #lhist is the newly proposed history samples
                                XT2 <: Union{Nothing,Vector{Int64}},
                                T, F, S}

    if typeof(lstem) == Nothing
        debugOnce("probSum for MhSampled1",l)
        a = probSum(fp,i,d,l,l,l,lstem,lhist)
        debugOnce("probSum for MhSampled2",a)
        a
    else
        debugOnce("probSum unVoid")
        probSum(fp,i,d,l,l,l,lstem,lhist)
    end
end

function  probSumWeight(fp::AbstractFinkel,
            ln::Int64) #location for neighborhood center
    1 #stub; in future, look at distribution of forward weights...
end

"""
    probSum for MhCompromise

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(
    p::FinkelParams{S,MhCompromise},
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::AbstractFinkel,
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Nothing,Int64},
                                XT2 <: Union{Nothing,Vector{Int64}},
                                T, F, S}

    myCenters = prodNeighborhood( l, d, fp.params.mh.rSub)
    result = 0.
    if typeof(lstem) == Nothing
        for c = myCenters
            result += probSum(fp,i,d,c,c,nothing,nothing,nothing) * probSumWeight(fp,c)
        end
    else

        for c = myCenters
            result += probSum(fp,i,d,c,c,l,lstem,lhist) * probSumWeight(fp,c)
        end
    end
    result
end

"""
    probSum for MhMultisampled

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(
    p::FinkelParams{S,MhMultisampled},
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::AbstractFinkel,
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Nothing,Int64},
                                XT2 <: Union{Nothing,Vector{Int64}},
                                T, F, S}

    myCenters = prodNeighborhood( l, d, fp.params.mh.rSub)
    result = 0.
    if typeof(lstem) == Nothing
        for c = myCenters
            result += probSum(fp,i,d,c,c,nothing,nothing,nothing) * probSumWeight(fp,c)
        end
    else

        for c = myCenters
            result += probSum(fp,i,d,c,c,l,lstem,lhist) * probSumWeight(fp,c)
        end
    end
    result
end


function clearOldProbs!(
    fp::AbstractFinkel,
    i::Int64,#current particle
    l::Int64,
    d::Int64)

    clearOldProbs!(fp,i,l,d,fp.params)#do nothing - mhsampled means it's just local
end

function clearOldProbs!(
    fp::AbstractFinkel,
    i::Int64,#current particle
    l::Int64,
    d::Int64,
    p::FinkelParams{S,MhSampled}) where {S}

    #do nothing - mhsampled means it's just local
end

function clearOldProbs!(
    fp::AbstractFinkel,
    i::Int64,#current particle
    l::Int64,
    d::Int64,
    p::FinkelParams{S,MhCompromise}) where {S}


    myCenters = prodNeighborhood( l, d, fp.params.mh.rSub)
    for λ = myCenters
        fp.prevProbs[l,i] = nothing
    end
end

function clearOldProbs!(
    fp::AbstractFinkel,
    i::Int64,#current particle
    l::Int64,
    d::Int64,
    p::FinkelParams{S,MhMultisampled}) where {S}


    myCenters = prodNeighborhood( l, fp, d, fp.params.mh.r)
    for λ = myCenters
        fp.prevProbs[l,i] = nothing
    end
end


function calcPrevProbs!(fp::FinkelParticles, d, n)

    for i in 1:n
        for l in 1:d
            fp.prevProbs[l,i] = Nullable{Float64}()#Nullable(probSum(fp,l,i,d))
        end
    end
end

function prodNeighborhood(l::Int64, #neighborhood over which we use historyTerms
                    d::Int64, #dimension/number of loci
                    neighbors=4)
    [mod1(i,d) for i in (l-div(neighbors,2)):(l-div(neighbors,2)+neighbors-1)]
    #[max(l-w,1):min(l+w,d)]
end


function mcmc!(fp::AbstractFinkel,i::Int64,steps::Int64)
    d = size(fp.base,1)
    for s in 1:steps
        order = randperm(d)
        for l in order
            p = sample(fp.ws[l])
            if fp.base[l,p] != fp.tip.particles[l,i]
                oldProbNull = fp.prevProbs[l,i]
                debugOnce("oldProbNull",oldProbNull)
                if oldProbNull == nothing
                    oldProb = probSum(i, d, l, fp)
                    debugOnce("oldProb wasnull",oldProb,i,l,s)
                    fp.prevProbs[l,i] = oldProb
                else
                    oldProb = oldProbNull
                    debugOnce("oldProb nonnull",oldProb,i,l,s)
                end
                newHistoryTerms = histTerms(l, p, fp)
                nProb = newProb = probSum(i, d, l, fp, p, newHistoryTerms)
                debugOnce("newProb ",oldProb,i,l,newProb,s)
                if fp.params.useForward != 0.
                    newProb = newProb / exp(fp.params.useForward * fp.logForwardDensities[l,p])
                    oldProb = oldProb / exp(fp.params.useForward * fp.logForwardDensities[l,fp.stem[l,i]])
                    debugOnce("useForward ",oldProb,i,l,newProb,s)
                end

                if (newProb < oldProb) && (newProb < rand() * oldProb)
                    fp.totalProb[l,i] += newProb / oldProb
                    debugOnce("reject totalProb ",oldProb,i,l,newProb,s)
                    continue #M-H rejection
                end
                #M-H accepted

                fp.numMhAccepts += 1
                fp.totalProb[l,i] += 1
                debugOnce("accept totalProb ",oldProb,i,l,newProb,fp.totalProb[l,i],s)
                fp.stem[l,i] = p
                fp.historyTerms[l,i,:] = newHistoryTerms

                fp.tip.particles[l,i] = fp.base[l,p]
                clearOldProbs!(fp,i,l,d)
                fp.prevProbs[l,i] = nProb
            end
        end
    end
end

function predictUpdate(prev::AbstractFinkel, y::Observation, nIter::Int64=15, debug=true)
    T = prev.params.algo
    if nIter>0
        fp = T(prev)
    else
        fp = T(prev,nothing) #nosteps - save time
    end
    reweight!(fp, y) #set ws
    replant!(fp) #set tip from base
    if nIter>0
        #Threads.@threads for i in 1:fp.tip.n #crashes... fix later
        @sync @distributed for i in 1:fp.tip.n #need to use SharedArrays... probably not too hard actually
        #for i in 1:fp.tip.n
            mcmc!(fp,i,nIter)
            if debug & ((i % 40)==0)
                print("Ran particle ", i, " ", nIter, "; mean tp = ", mean(fp.totalProb[:,1:i]), "\n")
            end
        end
    end
    fp
end

function resample(state::AbstractFinkel)
    state #do nothing
end
