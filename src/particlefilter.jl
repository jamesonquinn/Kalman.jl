#included in bkf.jl

abstract type AbstractParticleFilter end

mutable struct ParticleSet{T,F<:KalmanFilter} <: AbstractParticleFilter
    filter::F
    n::Int64
    particles::Array{T,2} #[location, particle] #or is it [particle, location]? Got error, confused.
    weights::ProbabilityWeights
end

Base.copy(pw::ProbabilityWeights) = deepcopy(pw) #Annoying that I have to do this explicitly because otherwise copy() isn't type-preserving

function Base.copy(ps::ParticleSet)
    ParticleSet(ps.filter,ps.n,copy(ps.particles),copy(ps.weights)) #Note that the filter is NOT copied but just referenced!
end

function Observation(ps::ParticleSet, n::Int64)
  Observation(ps.particles[:,n], ps.filter.z)
end

function Observation(p::Vector, om::LinearObservationModel)
  Observation(om.h * p + rand(MvNormal(om.r)))
end

mutable struct Resample
  s::Array{Int64}
end

function Resample(n::Int64)
  Resample(collect(1:n))
end

function length(r::Resample)
  length(r.s)
end

function ParticleSet(kf::KalmanFilter, n::Int64)
  d = toDistribution(kf)
  ParticleSet(kf, n, rand(d,n), ProbabilityWeights(ones(n),n))
end

#function ap(f::KalmanFilter,ps::Array{T,2})
  #don't force constructing a pset
#end
function getCurFuzzes(filt, prevParts, params)
  d, n = size(prevParts)
  hoodd = params.mh.r
  Σ = cov(prevParts;dims=2)
  replace_nan!(Σ)
  Σinv = Matrix{Float64}(undef,d,d)
  try
    Σinv = inv(Σ)
  catch
    Σ += Matrix(1e-6I,d,d)
    Σinv = inv(Σ)
    debug("no inverse but I fixed it")
  end
  μ = mean(prevParts;dims=2)
  localVariance = Array{Float64,2}(undef,n,d) #particle, locus
  basePrecision = Vector{Float64}(undef,d)
  diffs = prevParts .- μ
  for l in 1:d
    hood = prodNeighborhood(l,d,hoodd)
    hooddet = det(Σ[hood,hood])
    dethat = prod(Σ[i,i] for i in hood)
    basePrecision[l] = Σ[l,l] * (hooddet/dethat)^(1/hoodd)
    for p in 1:n
      dist = MvNormal(Matrix(1e-4I,hoodd,hoodd)) #hopefully this will be replaced
      try
        dist = MvNormal(Σ[hood,hood])
      catch
        d = length(hood)
        try
          dist = MvNormal(Σ[hood,hood] + Matrix(1e-4I,hood,hood))
          debug("Cholesky 3")
        catch
          debug("Cholesky 4")
        end
      end
      logRelDens = logpdf(dist,diffs[hood,p]) - logpdf(dist,zeros(hoodd))
      localVariance[p,l] = 1/(n-hoodd)/(1+exp(logRelDens))
    end
  end

  # [inv(Σinv + Diagonal([basePrecision[l] / params.overlap *
  #       prod(exp(localVariance[p,l]/params.mh.r^2)
  #       for λ in prodNeighborhood(l,d,params.mh.r))
  #     for l in 1:d]))
  #   for p in 1:n]

  result = Vector{Matrix{Float64}}(undef,n)
  for p in 1:n
    try
      result[p] = inv(Σinv + Diagonal([basePrecision[l] / params.overlap *
             prod(localVariance[p,λ]^(1/hoodd^2)
                 #localVar is inherently dimension hoodd, and we have hoodd copies, so divide by hoodd^2
             for λ in prodNeighborhood(l,d,hoodd))
           for l in 1:d]))
    catch
      debug("Problem in getCurFuzzes", p, basePrecision)
      result[p] = inv(Σinv + Diagonal([basePrecision[l] / params.overlap for l in 1:d]))
    end
  end
  result
end


function getCurFuzzQuick(filt, prevParts, params)
  d, n = size(prevParts)
  hoodd = params.mh.r
  Σ = cov(prevParts;dims=2)
  replace_nan!(Σ)
  Σinv = Matrix{Float64}(undef,d,d)
  try
    Σinv = inv(Σ)
  catch
    Σ += Matrix(1e-6I,d,d)
    Σinv = inv(Σ)
    debug("no inverse but I fixed it")
  end
  return(Σinv * ((n-hoodd)/ params.overlap)^(1/hoodd) )
end

function getNextFuzzes(filt, prevParts, params)
  d, n = size(prevParts)
  curFuzzes = getCurFuzzes(filt, prevParts, params)
  [propagateUncertainty(filt.f,prevParts[:,p:p],
      curFuzzes[p])
    for p in 1:n]
end
#
# function rejuvenateFuzz!(pset::ParticleSet) #not used, and BROKEN — newCenters should not be here
#   myparams = params(pset)
#   if myparams.rejuv != 0
#     prevParts = particleMatrix(pset)
#     d, n = size(prevParts)
#     filt = getbkf(pset)
#     newCtrs = newCenters(filt, prevParts)
#     fuzzes = getCurFuzzes(filt, prevParts, prev.params)
#
#     for j = 1:n
#       pset.particles[:,j] += newCtrs[:,j] + rand(
#               MvNormal(Matrix(Hermitian(fuzzes[j] * myparams.rejuv))))
#           #Matrix(Hermitian( :  ...Work, stupid!
#           #/4 : rejuv lightly, but pretend it's full. Not actually correct but meh.
#     end
#   end
# end
#
# function rejuvenateFullCloud(pset::ParticleSet,newPortion = .5) #not used, obsolete
#   L,M = size(pset.particles)
#   Σ = cov(pset.particles,pset.weights,2) #Technically, I should restrict this to banded part, but meh.
#   μ = vec(mean(pset.particles,pset.weights,2))
#   nOld = floor(Int, (1-newPortion) * M)
#   wmean = mean(pset.weights[1:nOld])
#   newParticles = typeof(pset.particles)(undef,L,M)
#   newParticles[:,1:nOld] = pset.particles[:,1:nOld]
#   newParticles[:,(nOld+1):M] = rand(MvNormal(μ,Σ),M-nOld)
#   ParticleSet(pset.filter,pset.n,newParticles,
#             ProbabilityWeights([pset.weights[1:nOld];[wmean for i in 1:(M-nOld)]]))
# end

function ap(pset::ParticleSet)
  ε=rand(noiseDistribution(pset.filter),pset.n)
  ParticleSet(pset.filter, pset.n, newCenters(pset.filter, pset.particles) + ε, pset.weights)
end


function ap(pset::ParticleSet, samp::Resample)
  n = length(samp)
  ε=rand(noiseDistribution(pset.filter),n)
  newParticleBases = Array{Float64}(undef,size(pset.particles,1), n)
  #println(size(ε))
  #println(size(newParticles))
  for i in 1:n

    newParticleBases[:,i] = newCenters(pset.filter, pset.particles[:,samp.s[i]])
  end
  ParticleSet(pset.filter, n, newParticleBases + ε, ProbabilityWeights(ones(n),n))
end

function reweight!(pset::ParticleSet, y::Observation)
  ws = pdf(obsNoiseDistribution(pset.filter),
                                pset.particles - pset.filter.z.h * repeat(y.y,outer=[1,size(pset.particles,2)])
                                )
  replace_nan!(ws)
  pset.weights = ProbabilityWeights(ws)
end

function mysample(::Val{:systematic}, set, w::ProbabilityWeights, n)
  if 0<sum(w)<Inf
    inc=sum(w)/n
    basecut = rand() * inc
    j = 0
    results = Vector{Float64}(undef,n)
    for (i, cursum) in enumerate(cumsum(w))
      while (basecut + inc * j) < cursum
        j += 1
        if j>n
          debugOnce("mysample overflow", i,j,cursum,n,sum(w), inc)
          break
        end
        results[j] = i
      end
    end
    for count in 1:length(results) #double check, yuck.
      while !(0<results[count]<=n)
        results[count] = sample(set,w)
      end
    end
  else #w is useless for some reason.
    results = sample(set,n)
    for count in 1:length(results) #double check, yuck.
      while !(0<results[count]<=n)
        results[count] = 1
      end
    end
  end
  results
end

function Resample(pset::ParticleSet)
  Resample(mysample(Val(:systematic),1:pset.n,pset.weights,pset.n)) #TODO: systematic (low-variance) resampling; Thrun/burgard/fox
end

mutable struct ParticleStep <: AbstractParticleFilter
    r::Resample
    p::ParticleSet
    y::Observation
end

function ParticleStep(pset::ParticleSet)
  o = Observation(Float64[])
  r = Resample(pset.n)
  ParticleStep(r, pset, o)
end

function ParticleStep(pset::ParticleSet, y::Observation, rejuv=.5)
  if rejuv>0
    pset = rejuvenate(pset, rejuv, true) #true for quickie
  end
  r = Resample(pset)
  p = ap(pset,r)
  reweight!(p,y)
  ParticleStep(r,p,y)
end

function ParticleStep(pstep::ParticleStep, y::Observation, rejuv=.5)
  ParticleStep(pstep.p,y, rejuv)
end

function predictUpdate(pstep::ParticleStep, y::Observation)
  ParticleStep(pstep, y)
end


function ParticleStep(kf::KalmanFilter, n::Int64)
  ParticleStep(ParticleSet(kf,n))
end

function particles(p::ParticleStep)
    p.particles
end
