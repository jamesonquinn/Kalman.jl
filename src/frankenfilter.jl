#included in bkf.jl

  abstract type AbstractParticleFilter end

type FrankenSet{T,F<:KalmanFilter} <: AbstractParticleFilter
    filter::F
    n::Int64
    hoodSize::Int64
    particles::Array{T,2} #[location, particle]
    weights::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}
end


  function toFrankenSet(kf::KalmanFilter, n::Int64, hoodSize)
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
    resampledParticles = resampleParticles(pset,samp,n)
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
    for i in 1:size(pset.weights,1)
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
end

  function FrankenStep(pset::FrankenSet)
    o = Observation([0.])
    r = FResample(pset.n, size(pset.weights,1))
    FrankenStep(r, pset, o)
  end

  function ParticleSet(pset::FrankenSet)
      r = FResample(pset)
      n = pset.n
      resampledParticles = resampleParticles(pset,r,n)
      ParticleSet(pset.filter, n, resampledParticles, ProbabilityWeights(ones(n)))
  end

  function doFrankenStep(pset::ParticleSet, fset::FrankenSet, y::Observation)
      p = ap(fset, fset.n, pset.particles)
      reweight!(p,y)
      FrankenStep(Vector{Resample}(),p,y)
  end

  function FrankenStep(pset::FrankenSet, y::Observation)
    r = FResample(pset)
    p = ap(pset,r)
    reweight!(p,y)
    FrankenStep(r,p,y) #NOTE: Resample happens at beginning, not end... so all my tests have been measuring wrong.
  end

#Two ways to progress:
#1: todayFrank = FrankenStep(yesterdayFrank, y)
#2: yesterdaySet = ParticleSet(yesterdayFrank.p); todayFrank = FrankenStep(yesterdaySet, yesterdayFrank, y)

#need to implement #2 as resample, predictupdate

function resample(unresampled::FrankenStep) #returns tuple (yesterdaySet, yesterdayFrank)
    (ParticleSet(unresampled.p), unresampled)
end

function predictupdate(stuff::Tuple, y::Observation)
    (yesterdaySet, yesterdayFrank) = stuff
    if isa(yesterdayFrank,FrankenStep)
        yesterdayFrank = yesterdayFrank.p
    end
    doFrankenStep(yesterdaySet, yesterdayFrank, y)
end

  function FrankenStep(pstep::FrankenStep, y::Observation)
    FrankenStep(pstep.p,y)
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
