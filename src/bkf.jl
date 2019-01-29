module bkf

using Statistics
using StatsBase
using Distributions
using Memoize
using Compat
using ForwardDiff
using LinearAlgebra
using DelimitedFiles
using Random

#for repeatableKF
using DataFrames, CSV, DataStructures, Dates
using NaNMath


if false
  using Distributed
  using Arrays #SharedArray
  if length(workers())<4
    addprocs(4)
  end
end

import Base.length

abstract type KalmanFilter end

abstract type LinearKalmanFilter <: KalmanFilter end

abstract type Model end

abstract type ObservationModel end

abstract type AbstractState end

mutable struct State{T} <: AbstractState
    x::Vector{T}
    p::Matrix{Float64}  #covariance
end

Base.:(==)(x1::State,x2::State) = x1.x==x2.x && x1.p == x2.p

mutable struct LinearModel <: Model
    a::Matrix #state transition matrix
    g::Matrix #process noise to state
    q::Matrix #process noise covariance (hopefully diagonal)
end
  function dimOf(lm::LinearModel)
      size(lm.a,1)
  end

## This is the first function a Kalman filter needs to
# implement: Apply its model to its state to generate a
# new state

function ap(f::LinearModel,x::State)
    x1 = f.a*x.x
    p1 = f.a*x.p*f.a' + f.g*f.q*f.g'
    State(x1,p1)
end

mutable struct Observation{T}
    y::Vector{T}
end

#  Base.convert(::Type{Observation},y) = Observation([y])

mutable struct LinearObservationModel <: ObservationModel
    h::Matrix #connects state to observations; hopefully eye
    r::Matrix #observation covariance; hopefully diagonal, but not necessarily eye (?)
end

function LinearObservationModel(h::Matrix)
    LinearObservationModel(Matrix(1.0I,size(h,1),size(h,1)),h)
end


mutable struct BasicKalmanFilter <: LinearKalmanFilter
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
function covs(kf::KalmanFilter,y::Observation)
    res = y.y - kf.z.h * kf.x.x
    ph = kf.x.p * kf.z.h'
    s = kf.z.h * ph + kf.z.r
    (res,ph,s)
end

function toDistribution(kf::BasicKalmanFilter)
  MvNormal(kf.x.x,((kf.x.p + kf.x.p') / 2))
end

function forwardDistribution(m::LinearModel,x::Vector,r::AbstractRange)
    #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
    #print("bbb",noiseMatrix(m)[r,r],"\n")
    #print("ccc",x[1:2],"\n")
    #print("ddd",x[r],"\n")

  MvNormal(x[r],noiseMatrix(m)[r,r])
end
function forwardDistribution(m::LinearModel,x::Vector,l::Int64)
    #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
    #print("bbb",noiseMatrix(m)[r,r],"\n")
    #print("ccc",x[1:2],"\n")
    #print("ddd",x[r],"\n")

  forwardDistribution(m,x,l:l)
end

function forwardDistribution(m::LinearModel,x::Float64,l::Int64)
    #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
    #print("bbb",noiseMatrix(m)[r,r],"\n")
    #print("ccc",x[1:2],"\n")
    #print("ddd",x[r],"\n")

  Normal(x,noiseMatrix(m)[l,l])
end

@memoize function noiseMatrix(m::LinearModel)
    m.g * m.q * m.g'
end

function noiseDistribution(kf::BasicKalmanFilter)
  MvNormal(noiseMatrix(kf.f))
end

function obsNoiseDistribution(kf::BasicKalmanFilter)
  MvNormal(kf.z.r)
end

function obsNoiseDistribution(kf::BasicKalmanFilter, fromI, toI)
  MvNormal(kf.z.r[fromI:toI,fromI:toI])
end

function newCenters(kf::BasicKalmanFilter, oldState)
    debug("old version nopers")
    kf.f.a * oldState
end


include("debug.jl")
include("lorenz.jl")
include("particlefilter.jl")
include("frankenfilter.jl")
include("finkel.jl")
include("fuzzfinkel.jl")

include("filter.jl")
include("unscented.jl")
include("kldivergence.jl")
include("repeatableKF.jl")
#include("kldivergencedtsgsdtg.jl")

end
