#included in bkf.jl

abstract type AbstractParticleFilter end

mutable struct FrankenSet{T,F<:KalmanFilter} <: AbstractParticleFilter
    filter::F
    n::Int64
    hoodSize::Int64
    particles::Array{T,2} #[location, particle]
    weights::Vector{ProbabilityWeights{Float64,Float64,Array{Float64,1}}}
end

#3
mutable struct FrankenStep <: AbstractParticleFilter
    r::Vector{Resample}
    p::FrankenSet
    y::Observation
end

function toFrankenSet(kf::KalmanFilter, n::Int64, hoodSize)
  d = toDistribution(kf)
  width = dimOf(kf.f)
  nHoods = div(width, hoodSize)
  FrankenSet(kf, n, hoodSize, rand(d,n),
                  [ProbabilityWeights(ones(n),Float64(n)) for i in 1:nHoods])
                  #maybe declare generic type for array to prevent over-specific?
                  #ProbabilityWeights[ProbabilityWeights(ones(n),Float64(n)) for i in 1:nHoods]) #?
end


function rejuvenate1(pset::FrankenSet,newPortion = .5)
  L,M = size(pset.particles)
  w = div(L, pset.hoodSize)
  newParticles = typeof(pset.particles)(undef,L,M)
  newWeights = typeof(pset.weights)(undef,w)

  nOld = floor(Int, (1-newPortion) * M)
  for nHood in 1:w
          hood = ((nHood-1)*pset.hoodSize+1):((nHood)*pset.hoodSize)
          Σ = cov(pset.particles[hood,:];dims=2, weights=pset.weights) #Technically, I should restrict this to banded part, but meh.
          μ = mean(pset.particles[hood,:];dims=2, weights=pset.weights)
          wmean = mean(pset.weights[nHood][1:nOld])
          newParticles[hood,1:nOld] = pset.particles[hood,1:nOld]
          newParticles[hood,(nOld+1):M] = rand(MvNormal(μ,Σ),M-nOld)
          newWeights[nHood] = ProbabilityWeights([pset.weights[nHood][1:nOld];
                                                  Vector(wmean,M-nOld)])
  end
  ParticleSet(pset.filter,pset.n,pset.hoodSize,
            newParticles,newWeights)
end

#function ap(f::KalmanFilter,ps::Array{T,2})
  #don't force constructing a pset
#end

function ap(pset::FrankenSet)
  ε=rand(noiseDistribution(pset.filter),pset.n)
  FrankenSet(pset.filter, pset.n, pset.hoodSize, newCenters(pset.filter, pset.particles) + ε, pset.weights)
end

#2.2
function ap(pset::FrankenSet, samp::Vector{Resample})
  n = pset.n
  w = length(samp)
  ε=rand(noiseDistribution(pset.filter),n)
  resampledParticles = Array{Float64}(undef,size(pset.particles,1), n)

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
              newCenters(pset.filter, resampledParticles) + ε,
              [ProbabilityWeights(ones(n),Float64(n)) for i in 1:size(pset.weights,1)])
end

#2.3
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
  s = size(pset.weights,1)
  result = Vector{Resample}(undef, s)
  for hood in 1:s
    try
      result[hood] = Resample(mysample(Val(:systematic),1:pset.n,pset.weights[hood],pset.n))#systematic (low-variance) resampling; Thrun/burgard/fox
    catch
      debug("mysample failure in FResample", hood, pset.n, sum(pset.weights[hood])) #TODO: actually fix this bug, don't just workaround
      result[hood] = Resample(sample(1:pset.n,pset.weights[hood],pset.n))
    end
  end
  result
end

#2.1
function FResample(n::Int64, hoods)
  [Resample(collect(1:n))
    for hood in 1:hoods]
end

Base.copy(f::FrankenStep) = FrankenStep(f.r, f.p, f.y)


function FrankenStep(pset::FrankenSet)
  o = Observation([0.])
  #print("\n In FrankenStep3 \n")
  r = FResample(pset.n, size(pset.weights,1))
  #print("\n In FrankenStep4 \n")
  res= FrankenStep(r, pset, o)
  #print("\n In FrankenStep5 \n")
  res
end

#2
function FrankenStep(pset::FrankenSet, y::Observation, rejuv=0.)
  if rejuv != 0
    pset = rejuvenate(pset)
  end
  r = FResample(pset)
  p = ap(pset,r)
  reweight!(p,y)
  res = FrankenStep(r,p,y) #NOTE: Resample happens at beginning, not end... so all my tests have been measuring wrong.
  res
end

#1
function FrankenStep(pstep::FrankenStep, y::Observation, rejuv=0.)
  FrankenStep(pstep.p,y, rejuv)
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

#0
function predictUpdate(f::Union{FrankenSet,FrankenStep}, y::Observation, rejuv::Float64=0.)
    FrankenStep(f,y,rejuv)
end
