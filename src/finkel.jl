#included in bkf.jl

include("delegatemacro.jl")

abstract AbstractSparseFilter <: KalmanFilter

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

abstract AbstractFinkel <: AbstractParticleFilter

type FinkelToe{T,F<:KalmanFilter} <: AbstractFinkel
    tip::ParticleSet{T,F}
end

function particleMatrix(fp::AbstractFinkel)
    fp.tip.particles
end

function getbkf(fp::AbstractFinkel)
    fp.tip.filter
end

type FinkelParticles{T,F<:KalmanFilter} <: AbstractFinkel
    tip::ParticleSet{T,F} #This is where we do MCMC and get the answer. Also holds the filter. It's got extra stuff; ignore.
    base::Array{T,2} #This is the raw 1-step progression from last time.
                #As with all similar arrays, size is [d,n]; that is, first index is space and second is particle.
        #[space, particle]
    prev::AbstractFinkel
    stem::Array{Int64,2} #tells which base each tip comes from
        #[space, particle]
    ws::Vector{WeightVec} #selection probabilities at each point; p(y_l|x^i_l)
        #[space][particle]
    #lps::Array{Nullable{Float64},2} #log probabilities from base to tip; size [n, n]
    means::Array{T,2} #base, progressed, without adding noise.
    prevProbs::Array{Float64,2} #sum over neighborhood history of local probability - avoid double calculation
    localDists::Array{Distribution,2} #for calculating probs
    totalProb::Vector{Float64} #convergence diagnostic
end

function FinkelParticles(prev::AbstractFinkel)
    d, n = size(particleMatrix(prev))
    print(d," dn ",n,"\n")
    tip = ap(prev.tip)
    stem = zeros(Int64,d,n)
    filt = getbkf(prev)
    for j = 1:n
        for i = 1:d
            stem[i,j] = j
        end
    end
    fp = FinkelParticles(tip,
                copy(tip.particles), #base
                prev,
                stem,
                Vector{WeightVec}(), #empty weights
                filt.f.a * particleMatrix(prev),
                Array{Float64,2}(d,n),
                Array{Distribution,2}(d,n),
                zeros(d),
                )
    getDists!(fp, d, n)
    calcPrevProbs!(fp, d, n)
    fp
end

function particles(fp::FinkelParticles)
    fp.tip.particles
end

function reweight!(fp::FinkelParticles, y::Observation)
    fp.ws = Vector{WeightVec}()
    diffs = fp.base - fp.tip.filter.z.h * repeat(y.y,outer=[1,fp.tip.n])  # fp.tip.f.z.h should probably be eye ?
    vars = diag(fp.tip.filter.z.r) #assumes fp.tip.f.z.r is diagonal
    for i = 1:size(fp.base,1)
        wvec = Vector{Float64}(fp.tip.n)
        for j = 1:fp.tip.n
            wvec[j] = exp(-diffs[i,j]^2 / vars[i] / 2)
        end
        push!(fp.ws,WeightVec(wvec))
    end
end

function getDists!(fp::FinkelParticles, d, n)
    for i in 1:n
        x = Vector(slice(fp.means,:,i))
        for l in 1:d
            fp.localDists[l,i] = forwardDistribution(getbkf(fp).f,x,max(l-1,1):min(l+1,d))
        end
    end
end


function probSum(fp::FinkelParticles,l::Integer,i::Integer,p::Integer,
            neighborhood::Range,
            centerOnly = false)
    ps = 0.
    prop = fp.tip.particles[neighborhood,i]
    hist = fp.stem[neighborhood,i]
    if p != fp.stem[l,i]
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
            centerOnly = false)
    probSum(fp,l,i,fp.stem[l,i],d,centerOnly)
end


function calcPrevProbs!(fp::FinkelParticles, d, n)

    print(d," ddnn ",n,"\n")




    for i in 1:n
        for l in 1:d
            fp.prevProbs[l,i] = probSum(fp,l,i,d)
        end
    end
end

function mcmc!(fp::FinkelParticles,i::Integer,steps::Integer)
    d = size(fp.base,1)
    for s in 1:steps
        order = randperm(d)
        for l in order
            p = sample(fp.ws[l])
            if p != fp.stem[l,i]
                neighorhood = max(l-1,1):min(l+1,d)
                oldSet = Set(fp.stem[i] for i in neighborhood)
                baseNewProb = newProb = probSum(fp,l,i,p,neighborhood)
                oldProb = fp.oldProbs[l,i]
                if !(p in oldSet)
                    newProb += probSum(fp,l,i,fp.stem[l,i],neighborhood,true)
                    oldProb += probSum(fp,l,i,p,neighborhood,true)
                end
                if newProb < oldProb && rand() > (newProb / oldProb)
                    fp.totalProb[l] += newProb / oldProb
                    continue #M-H rejection
                end
                #M-H accepted
                fp.totalProb[l] += 1
                fp.stem[l,i] = p
                fp.tip.particles[l,i] = fp.base[l,p]
                fp.prevProbs[l,i] = baseNewProb
            end
        end
    end
end

function FinkelParticles(prev::AbstractFinkel, y::Observation, nIter=5)
    fp = FinkelParticles(prev)
    reweight!(fp, y)
    for i in 1:fp.tip.n
        mcmc!(fp,i,nIter)
    end
    fp
end
