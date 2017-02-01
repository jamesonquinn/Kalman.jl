module bkf
  using Distributions
  using StatsBase
  import Base.length

  abstract KalmanFilter

  abstract LinearKalmanFilter <: KalmanFilter

  abstract Model

  abstract ObservationModel

  abstract AbstractState

  type State{T} <: AbstractState
      x::Vector{T}
      p::Matrix  #covariance
  end

  Base.(:(==))(x1::State,x2::State) = x1.x==x2.x && x1.p == x2.p

  type LinearModel <: Model
      a::Matrix #state transition matrix
      g::Matrix #process noise to state
      q::Matrix #process noise covariance (hopefully eye)
  end

  ## This is the first function a Kalman filter needs to
  # implement: Apply its model to its state to generate a
  # new state

  function ap(f::LinearModel,x::State)
      x1 = f.a*x.x
      p1 = f.a*x.p*f.a' + f.g*f.q*f.g'
      State(x1,p1)
  end

  type Observation{T}
      y::Vector{T}
  end

#  Base.convert(::Type{Observation},y) = Observation([y])

  type LinearObservationModel <: ObservationModel
      h::Matrix #connects state to observations
      r::Matrix #observation covariance
  end

  type BasicKalmanFilter <: LinearKalmanFilter
      x::State
      f::LinearModel
      z::LinearObservationModel
  end

  Base.copy(kf::KalmanFilter) = deepcopy(kf)

  ## This is the second of two functions a Kalman filter needs to implement
  # Return three matrices:
  # 1. The residual (y-Hx) (`res`)
  # 2. The cross-covariance between state and measurement (`ph`)
  # 3. The innovation covariance matrix
  function covs(kf::BasicKalmanFilter,y::Observation)
      res = y.y - kf.z.h * kf.x.x
      ph = kf.x.p * kf.z.h'
      s = kf.z.h * ph + kf.z.r
      (res,ph,s)
  end

  function toDistribution(kf::BasicKalmanFilter)
    MvNormal(kf.x.x,kf.x.p)
  end

  abstract AbstractParticle

  type Particle{T} <: AbstractParticle
      x::Vector{T}
      logWeight::Float64
  end

  function noiseDistribution(kf::BasicKalmanFilter)
    MvNormal(kf.f.g*kf.f.q*kf.f.g')
  end

  function obsNoiseDistribution(kf::BasicKalmanFilter)
    MvNormal(kf.z.r)
  end

  function Particle(x)
    Particle(x,1.)
  end

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

  function ap(pset::ParticleSet)
    noise=rand(noiseDistribution(pset.filter),pset.n)
    ParticleSet(pset.filter, pset.n, pset.filter.f.a * pset.particles, pset.weights)
  end

  function ap(pset::ParticleSet, samp::Resample)
    n = length(samp)
    ε=rand(noiseDistribution(pset.filter),n)
    newParticleBases = Array{Float64}(size(pset.particles)[1], n)
    #println(size(ε))
    #println(size(newParticles))
    for i in 1:n

      newParticleBases[:,i] = pset.filter.f.a * pset.particles[:,samp.s[i]]
    end
    ParticleSet(pset.filter, n, newParticleBases + ε, WeightVec(ones(n),n))
  end

  function reweight!(pset::ParticleSet, y::Observation)
    pset.weights = WeightVec(pdf(obsNoiseDistribution(pset.filter),
                                  pset.particles - repeat(y.y,outer=[1,size(pset.particles)[2]])
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


  include("filter.jl")
  include("unscented.jl")
end
