#included in bkf.jl

  abstract AbstractParticleFilter

type ParticleSet{T,F<:KalmanFilter} <: AbstractParticleFilter
    filter::F
    n::Int64
    particles::Array{T,2}
    weights::WeightVec
end

  function Observation(ps::ParticleSet, n::Int64)
    Observation(ps.particles[:,n], ps.filter.z)
  end

  function Observation(p::Vector, om::LinearObservationModel)
    Observation(om.h * p + rand(MvNormal(om.r)))
  end

type Resample
    s::Array{Int64}
end

  function Resample(n::Int64)
    Resample(collect(1:n))
  end

  function length(r::Resample)
    length(r.s)
  end

  function toParticleSet(kf::KalmanFilter, n::Int64)
    d = toDistribution(kf)
    ParticleSet(kf, n, rand(d,n), WeightVec(ones(n),n))
  end

  #function ap(f::KalmanFilter,ps::Array{T,2})
    #don't force constructing a pset
  #end

  function ap(pset::ParticleSet)
    ε=rand(noiseDistribution(pset.filter),pset.n)
    ParticleSet(pset.filter, pset.n, pset.filter.f.a * pset.particles + ε, pset.weights)
  end


  function ap(pset::ParticleSet, samp::Resample)
    n = length(samp)
    ε=rand(noiseDistribution(pset.filter),n)
    newParticleBases = Array{Float64}(size(pset.particles,1), n)
    #println(size(ε))
    #println(size(newParticles))
    for i in 1:n

      newParticleBases[:,i] = pset.filter.f.a * pset.particles[:,samp.s[i]]
    end
    ParticleSet(pset.filter, n, newParticleBases + ε, WeightVec(ones(n),n))
  end

  function reweight!(pset::ParticleSet, y::Observation)
    pset.weights = WeightVec(pdf(obsNoiseDistribution(pset.filter),
                                  pset.particles - pset.filter.z.h * repeat(y.y,outer=[1,size(pset.particles,2)])
                                  ))
  end

  function Resample(pset::ParticleSet)
    Resample([sample(pset.weights) for i in 1:pset.n])
  end

type ParticleStep <: AbstractParticleFilter
    r::Resample
    p::ParticleSet
    y::Observation
end

  function ParticleStep(pset::ParticleSet)
    o = Observation([0.])
    r = Resample(pset.n)
    ParticleStep(r, pset, o)
  end

  function ParticleStep(pset::ParticleSet, y::Observation)
    r = Resample(pset)
    p = ap(pset,r)
    reweight!(p,y)
    ParticleStep(r,p,y)
  end

  function ParticleStep(pstep::ParticleStep, y::Observation)
    ParticleStep(pstep.p,y)
  end

function particles(p::ParticleStep)
    p.particles
end
