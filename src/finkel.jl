#included in bkf.jl

include("delegatemacro.jl")

abstract type AbstractSparseFilter <: KalmanFilter end

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



abstract type FinkelParams{S<:SampleType} end

type UniformFinkelParams{S<:SampleType} <: FinkelParams{S}
    histPerLoc::Int64 #"H", in paper; num histories per locus to calculate likelihood product for
    radius::Int64 #neighborhood size. radius of 1 = magnitude 3.
    s::S #SampleType
end

function fparams(basedOn::Any)
    UniformFinkelParams(2,1,SampleUniform())
end

type FinkelParticles{T,F<:KalmanFilter,P<:FinkelParams} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    base::Array{T,2} #This is the raw 1-step progression from last time.
                #As with all similar arrays, size is [d,n]; that is, first index is space and second is particle.
        #[space, particle]
    prev::AbstractFinkel
    historyTerms::Array{Int64,3} #tells which history each particle came from
        #[space, sample, particle]
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
end

function fparams(fp::FinkelParticles)
    fp.params
end


function FinkelParticles(prev::AbstractFinkel,
                         )
    d, n = size(particleMatrix(prev))
    myparams = fparams(prev)
    h = myparams.histPerLoc
    print(d," dn ",n,"\n")
    base = ap(prev.tip).particles
    tipVals = copy(base)


    historyTerms = zeros(Int64,d,h,n)
    # for j = 1:n
    #     for i = 1:d
    #         for s = 1:h
    #             historyTerms[i,s,j] = j
    #     end
    # end

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
                Vector{ProbabilityWeights}(), #empty weights
                lps,
                histSampProbs,
                means, #means
                Array{Nullable{Float64},2}(d,n), #prevProbs
                localDists, #localDists
                zeros(d,n), #totalProb
                myparams
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
          fp.historyTerms[l,i] = p
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


function probSum(fp::FinkelParticles,
  l::Integer, #locus
  i::Integer, #current particle
  p::Integer, #history
            neighborhood::Range,
            centerOnly = false)
    ps = 0.
    prop = fp.tip.particles[neighborhood,i]
    hist = fp.historyTerms[neighborhood,i]
    if p != fp.historyTerms[l,i]
        if l==1
            place=1
        else
            place=2
        end
        hist[place] = p
        prop[place] = fp.base[l,p]
    end
    if centerOnly
        histSet = Set(p)
    else
        histSet = Set(hist)
    end
    for h in histSet
        aprob = pdf(fp.localDists[l,h],prop)
        for n in neighborhood
            if n != l
                aprob /= pdf(forwardDistribution(getbkf(fp).f,
                                        fp.means[n,h],
                                        n),
                                fp.tip.particles[n,i])
                        #prob, at point n, of going from particle h at time t-1 to
                        #value fp.tip.particles[n,i] at time t.
            end
        end
        ps += aprob
    end
    ps
end

function probSum(fp::FinkelParticles,l::Integer,i::Integer,p::Integer,d::Int64,
                centerOnly = false)

    probSum(fp,l,i,p,
            max(l-1,1):min(l+1,d),centerOnly
            )
end

function probSum(fp::FinkelParticles,l::Integer,i::Integer,d::Int64,
            centerOnly::Bool = false)
    probSum(fp,l,i,fp.historyTerms[l,i],d,centerOnly)
end


function calcPrevProbs!(fp::FinkelParticles, d, n)

    for i in 1:n
        for l in 1:d
            fp.prevProbs[l,i] = Nullable{Float64}()#Nullable(probSum(fp,l,i,d))
        end
    end
end

function mcmc!(fp::FinkelParticles,i::Integer,steps::Integer)
    d = size(fp.base,1)
    for s in 1:steps
        order = randperm(d)
        for l in order
            p = sample(fp.ws[l])
            if p != fp.historyTerms[l,i]
                neighborhood = max(l-1,1):min(l+1,d)
                oldSet = Set(fp.historyTerms[i] for i in neighborhood)
                baseNewProb = newProb = probSum(fp,l,i,p,neighborhood)
                oldProb = fp.prevProbs[l,i]
                if isnull(oldProb)
                    oldProb = probSum(fp, l, i, fp.historyTerms[l,i], neighborhood)
                    fp.prevProbs[l,i] = Nullable(oldProb)
                else
                    oldProb = get(oldProb)
                end
                if !(p in oldSet)
                    newProb += probSum(fp,l,i,fp.historyTerms[l,i],neighborhood,true)
                    oldProb += probSum(fp,l,i,p,neighborhood,true)
                end
                if newProb < oldProb && rand() > (newProb / oldProb)
                    fp.totalProb[l,i] += newProb / oldProb
                    continue #M-H rejection
                end
                #M-H accepted
                fp.totalProb[l,i] += 1
                fp.historyTerms[l,i] = p
                fp.tip.particles[l,i] = fp.base[l,p]
                fp.prevProbs[l,i] = baseNewProb
                for n in neighborhood
                    if n != l
                        fp.prevProbs[n,i] = Nullable{Float64}()
                    end
                end
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
