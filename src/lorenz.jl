abstract type AbstractLorenz <: Model end

mutable struct LorenzModel <: AbstractLorenz
  F::Float64 #forcing constant
  d::Int64 #Dimension
  s::Float64 #time step (substep)
  ss::Int64 #how many substeps per external step
  g::Matrix #process noise to state
  q::Matrix #process noise covariance (hopefully diagonal)
end

mutable struct EkfishLorenzModel <: AbstractLorenz
    F::Float64 #forcing constant
    d::Int64 #Dimension
    s::Float64 #time step
    volumeFactor::Float64 #shrinks propagated quasi-uncertainty; should result in particles neither crowded nor lonely over subspace of neighboring loci
end

mutable struct BasicLorenzFilter <: KalmanFilter
  x::State
  f::AbstractLorenz
  z::LinearObservationModel
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

  #p1 = J*diagm(0 =>diag(x.p))*J' + f.g*f.q*f.g' #"diagm(0 =>diag(" :don't allow off-diagonals to build up over multiple steps
  #ptrace = tr(p1)
  #maxScale = f.F^2* f.d
  #if ptrace > maxScale
      p1 = (f.F*3/8)^2 * Matrix(1.0I,f.d,f.d)
  #end
  #print("ap   ",p1[1:3,1:3],"\n"); print("\n","""print("ap   ",p1[1:3,1:3],"\n")""")
  State(x1,p1)
end



function propagateUncertainty(f::AbstractLorenz,x::Matrix,u) #u is a matrix of uncertainty
  #If ForwardDiff is working correctly, it magically does the looping over f.ss for us! Magic is good!
  J = ForwardDiff.jacobian(v -> newCenters(f,v), x)

  p1 = J*u*J' #+ f.g*f.q*f.g' #"diagm(0 =>diag(" :don't allow off-diagonals to build up over multiple steps
  p1
end

function propagateUncertainty(s) #s::ParticleSet
  propagateUncertainty(s.filter.f,s.particles,s.filter.z.r)
end



function toDistribution(lf::BasicLorenzFilter)
    debug("toDistribution")
    debug(typeof(lf),typeof(lf.x))
    pmat = lf.x.p
    mineig = minimum(eigvals((pmat + pmat') / 2))

    #print("toDistribution ",mineig,lf.x.p[1:3,1:3],"\n"); print("\n","""print("toDistribution ",mineig,lf.x.p[1:3,1:3],"\n")""")
    try
        MvNormal(lf.x.x,((lf.x.p + lf.x.p') / 2))
    catch
       mineig = minimum(eigvals((lf.x.p + lf.x.p') / 2))
       MvNormal(lf.x.x,((lf.x.p + lf.x.p') / 2) - (mineig*Matrix(1.0I,lf.f.d,lf.f.d)))
    end
end


function getNextFuzz(filt::BasicKalmanFilter, prevParts, r)
    nothing #Tells forwardDistribution to ignore this
end

function getNextFuzz(filt::BasicLorenzFilter, prevParts, r)
    center = mean(prevParts,dims=2)
    recentered = prevParts .- center
    correlations = recentered * recentered'
    d = size(correlations)[1]
    p = size(prevParts)[1]
    for i in 1:(d-r)
        for j in (i+r):d
            correlations[i,j] = 0.
            correlations[j,i] = 0.
        end
    end
    fuzzPerDim = (det(correlations)/p) ^ (1/d)
    fuzz = diagm(0 =>ones(d))

    #Make sure fuzz variance is not larger than full covariance in any main dimension
    factor = 1
    for i in 1:d #This is a quick-and-dirty algorithm; good enough I hope
        if fuzz[i,i] > correlations[i,i]
            factor = sqrt(fuzz[i,i]/correlations[i,i])
            for j in [mod1(i-1,d),mod1(i+1,d)]
                fuzz[j,j] *= factor
            end
            fuzz[i,i] = correlations[i,i]
        end
    end

    propagateUncertainty(filt.f,center,fuzz)
end


function forwardDistribution(m,x,fuzz::Nothing,r)
    forwardDistribution(m,x,r)
end

function forwardDistribution(m::AbstractLorenz,x::Vector,fuzz::Matrix,r::AbstractRange)
  #print("aaa",noiseMatrix(m)[1:2,1:2],"\n")
  #print("bbb",noiseMatrix(m)[r,r],"\n")
  #print("ccc",x[1:2],"\n")
  #print("ddd",x[r],"\n")
  try
    MvNormal(x[r],Matrix(Hermitian((noiseMatrix(m)+fuzz)[r,r])))
  catch y
    #
    print("XXXXX forwardDistribution error\n")
    print(y,"\n")
    print(Matrix(Hermitian((noiseMatrix(m)+fuzz)[r,r])),"\n")
    print("XXXXXXXX forwardDistribution\n")
    MvNormal(x[r],Matrix(1.0I,size(r)[1],size(r)[1]))
  end
end

function forwardDistribution(m::AbstractLorenz,x::Vector,r::AbstractRange)
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

function noiseDistribution(kf::BasicLorenzFilter)
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
    newState = oldState
    for n in 1:lm.ss
      oldState = copy(newState) #ideally, would be destructive each time except the first, to avoid reallocating memory
      for i in 1:lm.d
          newState[i,:] += lm.s * (
              (oldState[mod1(i+1, lm.d),:] - oldState[mod1(i-2, lm.d),:])
                  .* oldState[mod1(i-1, lm.d),:]
              - oldState[i,:] .+ lm.F
              )
          if i==lm.d && abs(newState[i,1]) > 100
              print("Error in newCenters: state too far from 0 (slipped out of attractor?)\n")
               # ,i,oldState,newState[i],"x\n",
               #                oldState[mod1(i+1, lm.d),:],"y\n",
               #                oldState[mod1(i-2, lm.d),:],"z\n",
               #                oldState[mod1(i-1, lm.d),:],"a\n",
               #                oldState[i,:],
               #                lm.s * (
               #                    (oldState[mod1(i+1, lm.d),:] - oldState[mod1(i-2, lm.d),:])
               #                        .* oldState[mod1(i-1, lm.d),:]
               #                    - oldState[i,:] + lm.F
               #                    ),
               #                "\n")
          end
      end
    end
    newState
end

function newCenters(lm::LorenzModel, oldState::Vector)
    newState = oldState
    for n in 1:lm.ss
      oldState = copy(newState)
      for i in 1:lm.d
        newState[i] += lm.s * (
            (oldState[mod1(i+1, lm.d)] - oldState[mod1(i-2, lm.d)])
                * oldState[mod1(i-1, lm.d)]
            - oldState[i] + lm.F
            )
      end
    end
    #print("newState", newState, "\n")
    newState
end
