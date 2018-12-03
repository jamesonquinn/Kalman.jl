abstract type AbstractLorenz <: Model end

  type LorenzModel <: AbstractLorenz
      F::Float64 #forcing constant
      d::Int64 #Dimension
      s::Float64 #time step
      g::Matrix #process noise to state
      q::Matrix #process noise covariance (hopefully diagonal)
  end

type EkfishLorenzModel <: AbstractLorenz
    F::Float64 #forcing constant
    d::Int64 #Dimension
    s::Float64 #time step
    volumeFactor::Float64 #shrinks propagated quasi-uncertainty; should result in particles neither crowded nor lonely over subspace of neighboring loci
end

  function dimOf(lm::AbstractLorenz)
      lm.d
  end

  ## This is the first function a Kalman filter needs to
  # implement: Apply its model to its state to generate a
  # new state

function ap(f::LorenzModel,x::State)
  x1 = newCenters(f,x.x)
  #J = ForwardDiff.jacobian(v -> newCenters(f,v), x.x)

  #p1 = J*diagm(diag(x.p))*J' + f.g*f.q*f.g' #"diagm(diag(" :don't allow off-diagonals to build up over multiple steps
  #ptrace = trace(p1)
  #maxScale = f.F^2* f.d
  #if ptrace > maxScale
      p1 = (f.F*3/8)^2 * eye(f.d)
  #end
  print("ap   ",p1[1:3,1:3],"\n")
  State(x1,p1)
end

function ap(f::EkfishLorenzModel,x::State)
  x1 = newCenters(f,x.x)
  J = ForwardDiff.jacobian(v -> newCenters(f,v), x.x)

  p1 = J*J' * x.volumeFactor
  State(x1,p1)
end


  function propagateUncertainty(f::AbstractLorenz,x::Matrix,u) #u is a matrix of uncertainty
      J = ForwardDiff.jacobian(v -> newCenters(f,v), x)

      p1 = J*u*J' #+ f.g*f.q*f.g' #"diagm(diag(" :don't allow off-diagonals to build up over multiple steps
      p1
  end

  function propagateUncertainty(s) #s::ParticleSet
      propagateUncertainty(s.filter.f,s.particles,s.filter.z.r)
  end

  type BasicLorenzFilter <: KalmanFilter
      x::State
      f::AbstractLorenz
      z::LinearObservationModel
  end


function toDistribution(kf::BasicLorenzFilter)
    mineig = minimum(eigvals((kf.x.p + kf.x.p') / 2))

    print("toDistribution ",mineig,kf.x.p[1:3,1:3],"\n")
    try
        MvNormal(kf.x.x,((kf.x.p + kf.x.p') / 2))
    catch
       mineig = minimum(eigvals((kf.x.p + kf.x.p') / 2))
       MvNormal(kf.x.x,((kf.x.p + kf.x.p') / 2) - (mineig*eye(kf.f.d)))
    end
end

function forwardDistribution(m,x,fuzz::Void,r)
    forwardDistribution(m,x,r)
end

function forwardDistribution(m::AbstractLorenz,x::Vector,fuzz::Matrix,r::Range)
  #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
  #print("bbb",noiseMatrix(m)[r,r],"\n")
  #print("ccc",x[1:2],"\n")
  #print("ddd",x[r],"\n")

  MvNormal(x[r],(noiseMatrix(m)+fuzz)[r,r])
end

function forwardDistribution(m::AbstractLorenz,x::Vector,r::Range)
  #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
  #print("bbb",noiseMatrix(m)[r,r],"\n")
  #print("ccc",x[1:2],"\n")
  #print("ddd",x[r],"\n")

  MvNormal(x[r],noiseMatrix(m)[r,r])
end

function forwardDistribution(m::AbstractLorenz,x::Vector,fuzz::Matrix,l::Int64)
  #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
  #print("bbb",noiseMatrix(m)[r,r],"\n")
  #print("ccc",x[1:2],"\n")
  #print("ddd",x[r],"\n")

  forwardDistribution(m,x,fuzz,l:l)
end

function forwardDistribution(m::AbstractLorenz,x::Vector,l::Int64)
  #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
  #print("bbb",noiseMatrix(m)[r,r],"\n")
  #print("ccc",x[1:2],"\n")
  #print("ddd",x[r],"\n")

  forwardDistribution(m,x,l:l)
end

function forwardDistribution(m::LorenzModel,x::Float64,fuzz::Matrix,l::Int64)
    #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
    #print("bbb",noiseMatrix(m)[r,r],"\n")
    #print("ccc",x[1:2],"\n")
    #print("ddd",x[r],"\n")

  Normal(x,(noiseMatrix(m)+fuzz)[l,l])
end
function forwardDistribution(m::LorenzModel,x::Float64,l::Int64)
    #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
    #print("bbb",noiseMatrix(m)[r,r],"\n")
    #print("ccc",x[1:2],"\n")
    #print("ddd",x[r],"\n")

  Normal(x,(noiseMatrix(m))[l,l])
end

@memoize function noiseMatrix(m::LorenzModel)
    m.g * m.q * m.g'
end

function noiseDistribution(kf::BasicLorenzFilter,m::LorenzModel)
  MvNormal(noiseMatrix(kf.f))
end

function noiseDistribution(kf::BasicLorenzFilter,m::LorenzModel)
  MvNormal(noiseMatrix(kf.f))
end

function noiseDistribution(kf::BasicLorenzFilter,m::LorenzModel)
  MvNormal(noiseMatrix(kf.f))
end

function obsNoiseDistribution(kf::BasicLorenzFilter)
  MvNormal(kf.z.r)
end

function obsNoiseDistribution(kf::BasicLorenzFilter, fromI, toI)
  MvNormal(kf.z.r[fromI:toI,fromI:toI])
end

function newCenters(kf::BasicLorenzFilter, oldState)
    newCenters(kf.f,oldState)
end

function newCenters(lm::LorenzModel, oldState)
    newState = copy(oldState)
    for i in 1:lm.d
        newState[i,:] += lm.s * (
            (oldState[mod1(i+1, lm.d),:] - oldState[mod1(i-2, lm.d),:])
                .* oldState[mod1(i-1, lm.d),:]
            - oldState[i,:] + lm.F
            )
        if abs(newState[i,1]) > 100
            print("too big ",i,oldState,newState[i],"x\n",
                            oldState[mod1(i+1, lm.d),:],"y\n",
                            oldState[mod1(i-2, lm.d),:],"z\n",
                            oldState[mod1(i-1, lm.d),:],"a\n",
                            oldState[i,:],
                            lm.s * (
                                (oldState[mod1(i+1, lm.d),:] - oldState[mod1(i-2, lm.d),:])
                                    .* oldState[mod1(i-1, lm.d),:]
                                - oldState[i,:] + lm.F
                                ),
                            "\n")
        end
    end
    newState
end

function newCenters(lm::LorenzModel, oldState::Vector)
    newState = copy(oldState)
    for i in 1:lm.d
        newState[i] += lm.s * (
            (oldState[mod1(i+1, lm.d)] - oldState[mod1(i-2, lm.d)])
                * oldState[mod1(i-1, lm.d)]
            - oldState[i] + lm.F
            )
    end
    #print("newState", newState, "\n")
    newState
end
