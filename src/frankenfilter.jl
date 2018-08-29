#included in bkf.jl

  abstract type AbstractParticleFilter end

type FrankenSet{T,F<:KalmanFilter} <: AbstractParticleFilter
    filter::F
    n::Int64
    hoodSize::Int64
    particles::Array{T,2} #[location, particle]
    weights::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}
end


  function FrankenSet(kf::KalmanFilter, n::Int64, hoodSize::Int64)
    d = toDistribution(kf)
    width = size(kf.f.a, 1)
    nHoods = div(width, hoodSize)
    FrankenSet(kf, n, hoodSize, rand(d,n),
                    [ProbabilityWeights(ones(n),Float64(n)) for i in 1:nHoods])
                    #maybe declare generic type for array to prevent over-specific?
                    #ProbabilityWeights[ProbabilityWeights(ones(n),Float64(n)) for i in 1:nHoods]) #?
  end

  #function ap(f::KalmanFilter,ps::Array{T,2})
    #don't force constructing a pset
  #end

  function ap(pset::FrankenSet)
    ε=rand(noiseDistribution(pset.filter),pset.n)
    FrankenSet(pset.filter, pset.n, pset.hoodSize, pset.filter.f.a * pset.particles + ε, pset.weights)
  end


function resampleParticles(pset::FrankenSet, samp::Vector{Resample}, n)
  FrankenSet(pset.filter,n,pset.hoodSize,
              rawResampleParticles(pset,samp,n),
              typeof(pset.weights)())
end

  function rawResampleParticles(pset::FrankenSet, samp::Vector{Resample}, n)
    resampledParticles = Array{Float64}(size(pset.particles,1), n)

    w = length(samp)
    for hood in 1:w
        for i in 1:n
            for j in 1:pset.hoodSize
                l = (hood-1)*pset.hoodSize + j
                resampledParticles[l,i] = pset.particles[l,samp[hood].s[i]]
            end
        end
    end
    resampledParticles
end

  function ap(pset::FrankenSet, samp::Vector{Resample})
    n = pset.n
    resampledParticles = rawResampleParticles(pset,samp,n)
    ap(pset, n, resampledParticles)
end

function ap(pset::FrankenSet, n::Int64, resampledParticles::Array{Float64})
    ε=rand(noiseDistribution(pset.filter),n)

    #newParticleBases =
    #println(size(ε))
    #println(size(newParticles))
    FrankenSet(pset.filter, n, pset.hoodSize,
                pset.filter.f.a * resampledParticles + ε,
                [ProbabilityWeights(ones(n),Float64(n)) for i in 1:size(pset.weights,1)])
  end

  function reweight!(pset::FrankenSet, y::Observation)
    nHoods = size(pset.weights,1)
    if nHoods == 0
        width = size(pset.filter.f.a, 1)
        nHoods = div(width,pset.hoodSize)
        pset.weights = Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}(nHoods)
    end

    for i in 1:nHoods
        toI = i*pset.hoodSize
        fromI = toI - pset.hoodSize + 1
        pset.weights[i] = ProbabilityWeights(
                  pdf(obsNoiseDistribution(pset.filter, fromI, toI),
                      pset.particles[fromI:toI,:] - pset.filter.z.h[fromI:toI,fromI:toI] * repeat(y.y[fromI:toI],outer=[1,size(pset.particles,2)])
                  ))
    end
  end

  function FResample(pset::FrankenSet)
    [Resample([sample(pset.weights[hood]) for i in 1:pset.n])
        for hood in 1:size(pset.weights,1)]
  end

    function FResample(n::Int64, hoods)
      [Resample(collect(1:n))
        for hood in 1:hoods]
    end

type FrankenStep <: AbstractParticleFilter
    r::Vector{Resample}
    p::FrankenSet
    y::Observation
    needsresample::Bool
end

function FrankenStep(r::Vector{Resample}, p::FrankenSet ,y::Observation)
    FrankenStep(r,p,y,true)
end

function FrankenStep(p::FrankenSet, needsresample::Bool) #dummy values for r and y
    FrankenStep(Vector{Resample}(),
                p,
                Observation(Float64[]),
                needsresample)
end

function FrankenStep(kf::KalmanFilter, n::Int64, hoodSize::Int64) #build initial state from model
    p = FrankenSet(kf,n,hoodSize)
    FrankenStep(p, false)
end

function FrankenStep(pset::FrankenSet) #do resample
    o = Observation(Float64[])
    r = FResample(pset.n, size(pset.weights,1))
    FrankenStep(r, pset, o)
end

function FrankenStep(pset::FrankenStep) #do resample if needed
    if pset.needsresample
        FrankenStep(pset.p) #do resample
    else
        FrankenStep(pset.p, false) #just build object with dummies
    end
end

function FrankenStep(pset::FrankenSet) #do resample
      r = FResample(pset)
      n = pset.n
      FrankenStep(resampleParticles(pset,r,n),false)
end

function resample(pset::FrankenSet) #alias
    FrankenStep(pset)
end

function FrankenStep(fset::FrankenStep, y::Observation) #ap and reweight
      p = ap(fset.p, fset.p.n, fset.particles)
      reweight!(p,y)
      FrankenStep(Vector{Resample}(),p,y)
end

function predictupdate(fset::FrankenStep, y::Observation) #alias
    FrankenStep(fset, y)
end

function FrankenStep(pset::FrankenSet, y::Observation) #call (resample if needed, ap, reweight)
    FrankenStep(pset, y, true)
end

Base.copy(pset::FrankenSet) = deepcopy(pset)

function FrankenStep(pset::FrankenSet, y::Observation, needsresample::Bool) #resample if needed, ap, reweight
    if needsresample
        r = FResample(pset)
        p = ap(pset,r)
    else
        r = Vector{Resample}()
        p = copy(pset)
    end
    reweight!(p,y)
    FrankenStep(r,p,y) #NOTE: Resample happens at beginning, not end... so all my tests have been measuring wrong.
end

#Two ways to progress:
#1: todayFrank = FrankenStep(yesterdayFrank, y)
#2: yesterdaySet = FrankenStep(yesterdayFrank); todayFrank = FrankenStep(yeasterdaySet, y)





function FrankenStep(pstep::FrankenStep, y::Observation)
    FrankenStep(pstep.p, y, pstep.needsresample)
end

function particles(p::FrankenStep)
    p.particles
end

function nlocs(f::FrankenSet)
    size(f.particles)[1]
end

function nlocs(f::FrankenStep)
    nlocs(f.p)
end
