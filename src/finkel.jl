#included in bkf.jl

include("delegatemacro.jl")

abstract type AbstractSparseFilter <: KalmanFilter end

DEFAULT_PRODUCT_RADIUS = 1
DEFAULT_RADIUS_FRINGE = 0
DEFAULT_HISTPERLOC = 7

type SparseKF <: AbstractSparseFilter
    f::BasicKalmanFilter
    neighbors::Vector{Vector{Int16}}
end

canDebug = false
if canDebug

    verbose = false
    function debug(a...)
        global verbose
        if verbose
            for i=a
                print(i," ")
            end
            print("\n")
        end
    end
    function setVerbose(v=true)
        global verbose = v
    end
    debugDone = Dict{String,Bool}()
    function debugOnce(a...)
        global debugDone
        if !get(debugDone,a[1],false)
            debugDone[a[1]] = true
            for i=a
                print(i," ")
            end
            print("\n")
        end
    end

    function resetOnce()
        global debugDone = Dict{String,Bool}()
    end


else

    #stubs
    function debug(a...)
    end
    function debugOnce(a...)
    end
    function resetOnce()
    end

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

abstract type SampleType end #how to sample particle histories

type SampleUniform <: SampleType
end

type SampleLog <: SampleType
    inflectionPoint::Float64 #Distance below max log prob where probability slope bends
    factor::Float64 #how much probability slope bends
end

function SampleLog()
    SampleLog(5., 5.) #defaluts
end

type SampleLogM <: SampleType
    inflectionPoint::Float64 #Distance below max log prob where probability slope bends
    factor::Float64 #how much probability slope bends
    minSpan::Float64 #minimum distance from max to min
end

function SampleLogM(inflectionPoint::Float64,factor::Float64)
    SampleLog(inflectionPoint,factor,inflectionPoint*1.5)
end

abstract type MhType end #how to calculate MH acceptance probability

type MhLocal <: MhType #Restricted product, unrestricted sum
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
end

type MhSampled <: MhType #Neighborhood product, sum using locally sampled history
    r::Int64 #neighborhood size for product. radius of 1 = magnitude 3.
    histPerLoc::Int64 #history samples
end

function MhSampled(histPerLoc)
    MhSampled(DEFAULT_PRODUCT_RADIUS, histPerLoc)
end

type MhMultilocus <: MhType #Variable neighborhood product, sum using histories sampled for each center
    neighbors::Int64 #neighborhood size for product.
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

function MhCompromise(histPerLoc) #actually, total histories per attempt, rounded to nearest
    MhCompromise(DEFAULT_PRODUCT_RADIUS, round(Int64,histPerLoc/
                        (DEFAULT_PRODUCT_RADIUS - DEFAULT_RADIUS_FRINGE)),
                DEFAULT_PRODUCT_RADIUS - DEFAULT_RADIUS_FRINGE)
end

type FinkelParams{S<:SampleType,MH<:MhType}
    s::S
    mh::MH
    useForward::Float64
end


abstract type AbstractFinkel <: AbstractParticleFilter end

type FinkelToe{T,F<:KalmanFilter} <: AbstractFinkel
    tip::ParticleSet{T,F}
    params::FinkelParams
end

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
            useForward::Float64=1.)
    FinkelParams(SampleUniform(),MhSampled(radius, histPerLoc),
                useForward)
end

#uniform, compromise
function fparams(
            radius::Int64 ,
            histPerLoc::Int64 ,
            rSub::Int64,
            useForward::Float64=1.
            )
    FinkelParams(SampleUniform(),
                MhCompromise(radius, histPerLoc, rSub),
                useForward)
end

#log, sampled
function fparams(histPerLoc::Int64,
            inflectionPoint::Float64,
            factor::Float64,
            useForward::Float64=1.)
    FinkelParams(SampleLog(inflectionPoint,factor),MhSampled(DEFAULT_PRODUCT_RADIUS, histPerLoc),
                useForward)
end

#log, compromise
function fparams(inflectionPoint::Float64,
            factor::Float64,

            radius::Int64 ,
            histPerLoc::Int64 ,
            rSub::Int64,
            useForward::Float64=1.
            )
    FinkelParams(SampleLog(inflectionPoint,factor),
                MhCompromise(radius, histPerLoc, rSub),
                useForward)
end



#uniform, sampled is default
function fparams(basedOn::Any)
    fparams()
end

type FinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
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
        #[space, particle]
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

function fparams(fp::AbstractFinkel)
    fp.params
end

function FinkelParticles(prev::AbstractFinkel)
    myparams = fparams(prev)
    FinkelParticles(prev,myparams)
end

function FinkelParticles(prev::AbstractFinkel, nosteps::Void)
    myparams = fparams(prev)
    FinkelParticles(prev,myparams,nosteps)
end

function FinkelParticles(prev::AbstractFinkel,
                         myparams::FinkelParams)
    d, n = size(particleMatrix(prev))
    h = myparams.mh.histPerLoc
    base = ap(prev.tip).particles
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

    filt = getbkf(prev)
    means = newCenters(filt, particleMatrix(prev))
    histSampProbs = Array{ProbabilityWeights,2}(d,n)
    localDists = Array{Distribution,2}(d,n)
    logForwardDensities = Array{Float64,2}(d,n)
    for ph = 1:n
        for l = 1:d #Ideally we could somehow call forwardDistribution just once but meh
            localDists[l,ph] = forwardDistribution(filt.f,
                                    means[:,ph],
                                    l)
        end
    end

    lps = Array{Float64,3}(d,n,n)
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
    getDists!(fp, d, n)
    #calcPrevProbs!(fp, d, n)
    fp
end


function FinkelParticles(prev::AbstractFinkel,
                         myparams::FinkelParams,
                         nosteps::Void)
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

function getSampProbs(logpdfs,
                    pf::Int64,
                    s::SampleUniform)
    ProbabilityWeights(ones(logpdfs))
end

function getSampProbs(lpdfs,
                    pf::Int64,
                    s::SampleLog)
    logpdfs = copy(lpdfs)
    mn, mx = extrema(logpdfs)
    mn -= 1
    #logpdfs[pf] = mn #avoid (over-juicing low-density futures by choosing the past they're conditional on)
    if (mx - mn) < s.inflectionPoint
        ProbabilityWeights(logpdfs - mn)
    else
        v = logpdfs - mn
        inflec = mx - mn - s.inflectionPoint
        for i = 1:length(v)
            if v[i] > inflec
                v[i] += (v[i] - inflec) * s.factor
            end
        end
        ProbabilityWeights(v)
    end
end

function particles(fp::FinkelParticles)
    fp.tip.particles
end

function reweight!(fp::FinkelParticles, y::Observation)
    d = size(fp.base,1)
    fp.ws = Vector{ProbabilityWeights}(d)
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
        wvec = exp.(-diffs[l,:].^2 / vars[l] / 2 ) / forwardProb # ./ if forwardProb is vector
        fp.ws[l] = ProbabilityWeights(wvec)
    end
end

function replant!(fp::FinkelParticles)
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
        neighborhood::Vector,
        l::Void,
        lstem::Void)
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

function getSampProb(fp::FinkelParticles,
        i::Int64,#current particle
        l::Int64, #location for neighborhood center
        lstem::Int64)

    fp.histSampProbs[l,i][lstem] / sum(fp.histSampProbs[l,i])
end

function getSampProb(fp::FinkelParticles,
        i::Int64,#current particle
        ln::Int64, #location for neighborhood center
        lr::Void = nothing, #location to replace
        lstem::Void = nothing)
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

function getSampProb(fp::FinkelParticles     ,#{T,F,FinkelParams{S,MhMultisampled}},
        i::Int64,#current particle
        ln::Int64, #location for neighborhood center
        myneighborhood::Vector,## was UnitRange{Int64},
        lr::Union{Int64,Void}, #location to replace
        lstem::Union{Int64,Void}) #where {T,F,S}

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
function probSum(fp::FinkelParticles{},
    i::Int64,#current particle
    d::Int64,
    locn::Int64, #location for neighborhood center
    lh::Int64, #location for history samples
    lr::T = nothing, #location to replace
    lstem::T = nothing,
    lhist::Void = nothing) where T <: Union{Void,Int64}

    myneighborhood = prodNeighborhood( locn, fp, d, fp.params.mh.r)
    prob = 0.
    for h in fp.historyTerms[lh,i,:]
        prob += probSum(fp, i, h, myneighborhood,
            lr, lstem
            ) / getSampProb(fp, i, locn, myneighborhood, lr, lstem)
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
        lhist = Vector(fp.historyTerms[lh,i,:])
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
    lhist::XT2 = nothing) where {XT <: Union{Void,Int64}, #lhist is the newly proposed history samples
                                XT2 <: Union{Void,Vector{Int64}},
                                T, F, S}

    if typeof(lstem) == Void
        debugOnce("probSum for MhSampled1",l)
        a = probSum(fp,i,d,l,l,lstem,lstem,lhist)
        debugOnce("probSum for MhSampled2",a)
        a
    else
        debugOnce("probSum unVoid")
        probSum(fp,i,d,l,l,l,lstem,lhist)
    end
end

function  probSumWeight(fp::FinkelParticles{},
            ln::Int64) #location for neighborhood center
    1 #stub; in future, look at distribution of forward weights...
end

"""
    probSum for MhCompromise

    multiply likelihood from history samples for locus lh to neighborhood of particle i.
    If 0<l<=d, use lstem instead of current value at location l.
"""
function probSum(
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::FinkelParticles{T,F,FinkelParams{S,MhCompromise}},
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Void,Int64},
                                XT2 <: Union{Void,Vector{Int64}},
                                T, F, S}

    myCenters = prodNeighborhood( l, fp, d, fp.params.mh.rSub)
    result = 0.
    if typeof(lstem) == Void
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
    i::Int64,#current particle
    d::Int64,
    l::Int64,
    fp::FinkelParticles{T,F,FinkelParams{S,MhMultisampled}},
    lstem::XT = nothing,
    lhist::XT2 = nothing) where {XT <: Union{Void,Int64},
                                XT2 <: Union{Void,Vector{Int64}},
                                T, F, S}

    myCenters = prodNeighborhood( l, fp, d, fp.params.mh.rSub)
    result = 0.
    if typeof(lstem) == Void
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
    fp::FinkelParticles{T,F,FinkelParams{S,MhSampled}},
    i::Int64,#current particle
    l::Int64,
    d::Int64) where {T, F, S}

    #do nothing - mhsampled means it's just local
end

function clearOldProbs!(
    fp::FinkelParticles{T,F,FinkelParams{S,MhCompromise}},
    i::Int64,#current particle
    l::Int64,
    d::Int64) where {T, F, S}


    myCenters = prodNeighborhood( l, fp, d, fp.params.mh.rSub)
    for λ = myCenters
        fp.prevProbs[l,i] = nothing
    end
end

function clearOldProbs!(
    fp::FinkelParticles{T,F,FinkelParams{S,MhMultisampled}},
    i::Int64,#current particle
    l::Int64,
    d::Int64) where {T, F, S}


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
                    neighbors=4)
    [mod1(i,d) for i in (l-div(neighbors,2)):(l-div(neighbors,2)+neighbors-1)]
    #[max(l-w,1):min(l+w,d)]
end


function mcmc!(fp::FinkelParticles,i::Int64,steps::Int64)
    d = size(fp.base,1)
    for s in 1:steps
        order = randperm(d)
        for l in order
            p = sample(fp.ws[l])
            if fp.base[l,p] != fp.tip.particles[l,i]
                oldProbNull = fp.prevProbs[l,i]
                debugOnce("oldProbNull",oldProbNull)
                if isnull(oldProbNull)
                    oldProb = probSum(i, d, l, fp)
                    debugOnce("oldProb wasnull",oldProb,i,l,s)
                    fp.prevProbs[l,i] = Nullable(oldProb)
                else
                    oldProb = get(oldProbNull)
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

function FinkelParticles(prev::AbstractFinkel, y::Observation, nIter=15, debug=true)
    if nIter>0
        fp = FinkelParticles(prev)
    else
        fp = FinkelParticles(prev,nothing) #nosteps - save time
    end
    reweight!(fp, y) #set ws
    replant!(fp) #set tip from base
    if nIter>0
        for i in 1:fp.tip.n
            mcmc!(fp,i,nIter)
            if debug & ((i % 40)==0)
                print("Ran particle ", i, "; mean tp = ", mean(fp.totalProb[:,1:i]), "\n")
            end
        end
    end
    fp
end

function predictupdate(prev::AbstractFinkel, y::Observation, nIter::Int64)
    FinkelParticles(prev,y,nIter)
end

function resample(state::AbstractFinkel)
    state #do nothing
end
