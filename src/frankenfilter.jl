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


  function ap(pset::FrankenSet, samp::Vector{Resample})
    n = pset.n
    w = length(samp)
    ε=rand(noiseDistribution(pset.filter),n)
    resampledParticles = Array{Float64}(size(pset.particles,1), n)

    for hood in 1:w
        for i in 1:n
            for j in 1:pset.hoodSize
                l = (hood-1)*pset.hoodSize + j
                resampledParticles[l,i] = pset.particles[l,samp[hood].s[i]]
            end
        end
    end

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

  function FrankenStep(pset::FrankenSet, y::Observation)
    r = FResample(pset)
    p = ap(pset,r)
    reweight!(p,y)
    FrankenStep(r,p,y) #NOTE: Resample happens at beginning, not end... so all my tests have been measuring wrong.
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

function predictupdate(f::Union{FrankenSet,FrankenStep}, y::Observation)
    FrankenStep(f,y)
end
