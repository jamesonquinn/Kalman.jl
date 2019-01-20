function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

cd(joinpath(splitdir(@__FILE__)[1],"../output"))

ppath("/Users/chema/Dropbox/Kalman.jl/src")
ppath("/Users/chema/Dropbox/")
ppath("/Users/chema/mydev/Gadfly.jl/src")

# macro load(pkg)
#   quote
#     if isdefined(Symbol($(string(pkg))))
#       #try
#       #  reload(ucfirst(string(pkg)))
#       #catch
#         reload($(string(pkg)))
#       #end
#     else
#       import $pkg
#     end
#   end
# end
#
#
# @load bkf

using Revise
using Distributions
using bkf
using DelimitedFiles
using Random
using LinearAlgebra
using PlotlyJS

myWritecsv(io, a; opts...) = writedlm(io, a, ','; opts...)

#room for atom error display








valfname = "lr96fixedest10.2.csv"

function trsp(v)
    reshape(v,(1,:))
end
fname = "lorenz_graphjunk.csv"
#@load bkf

using bkf


NEIGHBORHOOD_SIZE = 4
IDEAL_SAMPLES = 1000

histPerLocs = [45,30,9]
nIters = [80,0,160,20,40]
useForwards = [1.,.5,0.]
#
if true #false for quickie test
    nParticles = [ #d,nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot,
                    #max sampType/mhType, max useForward
              #(60,80,  80^2   *10,div(80^2*2, 1), 4,2,20,1),
              #
              (20,200,40^2      ,div(200^2,5),1,1,20,1,1,1),
              ]
else
# nParticles = [(5,25,5,40,5,10,1), #nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot
#             (100, 100,20,10,5,10,3),
#             (500,250,10,7,3,5,2),
#             (1000,10,10,4,3,3,2)]
    nParticles = [ #d,nfp,npf,nfapf,reps,max nIters slot, steps,max histPerLoc slot,
                    #max sampType/mhType, max useForward
              #(60,80,  80^2   *10,div(80^2*2, 1), 4,2,20,1),
              #
              (9,80,4^2      ,div(80^2,5),1,1,20,1,1,1),
              ]
end

#Small numbers for quicker test
# nParticles = [(5,25,5,5,5,4), #nfp,npf,nfapf,reps,max nIters slot, steps
#             (100, 10000,2000,5,5,4)]
global sampTypes = [bkf.SampleLog(5.,5.), bkf.SampleUniform()]
global mhTypes = [bkf.MhSampled, bkf.MhCompromise]
global sampTypes = [bkf.SampleUniform(), bkf.SampleLog(5.,5.)]
global mhTypes = [bkf.MhSampled]
#mhTypes = []
global reps = max([_np[4] for _np in nParticles]...)#max of reps above
lnIters = length(nIters)
# finkelmeand = zeros(lnIters,ln,reps)
# frankenmeand = zeros(lnIters,ln,reps)
# partmeand = zeros(lnIters,ln,reps)
# idealmeand = zeros(lnIters,ln,reps)

global lnParts = length(nParticles)
global finkelkl = zeros(lnParts,1,reps)
global frankenkl = zeros(lnParts,1,reps)
global partkl = zeros(lnParts,1,reps)
global idealkl = zeros(lnParts,1,reps)
global np = 1
global nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc = nParticles[np]
global width = 1
try
    global mhType = mhTypes[1]
catch
    global mhTypes = bkf.MhSampled
end
global sampType = sampTypes[1]







global forcingF = 8.
#d is set by nParticles above
global timeSuperStep = .1
global numSteps = 15
global timeStep = timeSuperStep/numSteps
if timeStep > .01
  Val("Error here! increase numSteps or decrease timeSuperStep.")
end
global processNoiseVar = 0.0001 #Is this good? Needs testing.
global measurementNoiseVar = 0.1 #Again, ???
global useMeasurementNoiseVar = false #So ignore above line.
global initialvar = 0.4 #leaves room for early progress



global basenoise = .04
global highnoise = 4
global highgap = 5

global overlapFactor = 2.



for _np in 1:1 #jfc this is stupid but trust me don't change it

    global d,nfp,npf,nfapf,reps,lnIters,T,lnHistPerLoc,lnParams,lnForward = nParticles[_np]


    x0 = bkf.State(ones(d)*forcingF,initialvar*Matrix(1.0I,d,d))

    g = Matrix(1.0I,d,d) * processNoiseVar

    q = Matrix(1.0I,d,d)

    f = bkf.LorenzModel(forcingF,d,timeStep,numSteps,g,q)#,overlapFactor)

    if useMeasurementNoiseVar
        #Uniform noise
        r = Matrix(1.0I,d,d)*measurementNoiseVar
    else
        #low noise, with exceptions
        noisevec = fill(basenoise,d)
        noisevec[1:highgap:d] .= highnoise
        r = Diagonal(noisevec)
    end

    #var35 = bkf.meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = bkf.LinearObservationModel(Array(r))
    var35 = bkf.meanvarlocs(zeros(d), inv(r), 3:5)[2]

    kf0 = bkf.BasicLorenzFilter(x0,f,z)

    r = 1


    fpf = bkf.ParticleSet(kf0,nfp)
    pf = bkf.ParticleSet(kf0,npf)
    global fapf = bkf.toFrankenSet(kf0,nfapf,NEIGHBORHOOD_SIZE)

    pf1 = bkf.ParticleSet(kf0,1)

    t = collect(0:T)


    global truth = Vector{bkf.ParticleSet}(undef,0)#length(t))
    global observations = Vector{bkf.Observation}(undef,0)#length(t))
    global kfs = Vector{bkf.BasicLorenzFilter}(undef,0)#length(t))
    global pfs = Vector{bkf.ParticleStep}(undef,0)#length(t))
    global faps = Vector{bkf.FrankenStep}(undef,0)
    global fips = Vector{bkf.AbstractFinkel}(undef,0)

    global kf = kf0
    push!(kfs, kf)
    global ps = bkf.ParticleStep(pf)
    push!(pfs, ps)
    global fap = bkf.FrankenStep(fapf)
    #print(fap.p.particles[1:2,1:2], " fap\n")
    #print(bkf.musig(fap)[2][1:2,1:2], " fapμΣ\n")
    global fapSamps = copy(fap) #null resample
    push!(faps, fap)

    push!(truth, pf1)
    push!(observations, bkf.Observation(pf1,1))

    global nIter = 0
    #####
    global finalDist = bkf.toDistribution(kfs[end])
    pfinal = params(finalDist)
    # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
    # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
    # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
    # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

    #print("\nideal:"); print("\n","""print("\nideal:")""")
    global rsamps = rand(finalDist,IDEAL_SAMPLES)

    global kls = bkf.kl2(finalDist,rsamps)
    global sqe = bkf.sqerr(truth[end].particles[:,1],rsamps)
    global mvl = bkf.meanvarlocs(kfs[end],3:5)
    global idealkl[np,width,r] = kls[1]







    global r = 1

    for _r in 1:reps
        global r = _r



        #print("\nb15"); print("\n","""print("\nb15")""")
    ######

        i = 2

        for i in 2:length(t)-1

            open( fname,  "a") do outfile
                push!(truth, bkf.ap(truth[i-1]))
                push!(observations, bkf.Observation(truth[i],1))
                #print("\nideal:"); print("\n","""print("\nideal:")""")
                kf2 = bkf.predict(kf)
                kf = bkf.update(kf2,observations[i])
                push!(kfs, kf)

                print("\npart:")
                parttime = (@timed ps = bkf.ParticleStep(ps,observations[i]))[2]
                push!(pfs, ps)
                print("\nfranken:")
                franktime = (@timed begin

                    fap = bkf.FrankenStep(fapSamps, observations[i])

                    fapSamps = copy(fap)
                end)[2]

                push!(faps, fap)

                dif = kfs[end].x.p - kfs[end].x.p'
                mean(dif)
                mean(kfs[end].x.p)
                global finalDist = bkf.toDistribution(kfs[end])
                # finkelmeand[width,r] = mean(log.(bkf.pdf(finalDist,fps[end].tip.particles)))/d
                # partmeand[width,r] = mean(log.(bkf.pdf(finalDist,pfs[end].p.particles)))/d
                # frankenmeand[width,r] = mean(log.(bkf.pdf(finalDist,faps[end].p.particles)))/d
                # idealmeand[width,r] = mean(log.(bkf.pdf(finalDist,rand(finalDist,50))))/d

                rsamps = rand(finalDist,IDEAL_SAMPLES)
                kls = bkf.kl2(finalDist,rsamps)
                sqe = bkf.sqerr(truth[i].particles[:,1],rsamps)
                mvl = bkf.meanvarlocs(kfs[end],3:5)
                idealkl[np,width,r] = kls[1]
            end
            #print("\nb")
            print("Finkel ",r, " ... ",nfp); print("\n","""print("Finkel ",r, " ... ",nfp)""")
            ni, nhist, nIter, i = 1,1,1,2
        end
        global useforward, mhType, sampType, ni = [useForwards[1], mhTypes[1], sampTypes[1], 1]
        1
        for useForward in useForwards[1:lnForward]
            for mhType in mhTypes[1:min(length(mhTypes),lnParams)]
                for sampType in sampTypes[1:min(length(sampTypes),lnParams)]
                    for ni = 1:lnIters
                        nIter = nIters[ni]
                        print("nIter ",nIter)
                        for nhist = 1:lnHistPerLoc
                            histPerLoc = histPerLocs[nhist]


                            global fp = bkf.FinkelToe(fpf,bkf.FinkelParams(sampType,mhType(histPerLoc),useForward,
                                  overlapFactor,bkf.FuzzFinkelParticles,0.25)) #fparams(histPerLoc,2))
                            global ft = fp
                            #sampType = "sampled..uniform"
                            #fps = Vector{bkf.AbstractFinkel}(0)#length(t))
                            #push!(fps, fp)
                            print("\nFinkel:",mhType," ",sampType," ",np," ",width," ",r," ",nIter," ",i," ",histPerLoc," \n")
                            kls = bkf.kl2(finalDist,fp.tip.particles)
                            sqe = bkf.sqerr(truth[i].particles[:,1],fp.tip.particles)
                            finkelkl[np,width,r]  = kls[1]
                            print("\n  hm1:")
                            for i in 2:length(t)-1

                                open( fname,  "a") do outfile
                                    #
                                    finalDist = bkf.toDistribution(kfs[i])
                                    if true #false to time
                                      fp = bkf.predictUpdate(fp, observations[i], nIter)
                                      finktime = 0.
                                    else
                                      finktime = (@timed fp = bkf.predictUpdate(fp, observations[i], nIter))[2]
                                    end
                                    push!(fips,fp)
                                end
                            end
                        end

                    end #Mh
                end #Samp
            end
        end


    end
end




[idealkl,finkelkl,frankenkl,partkl

]

function evograph()
    #
    if true #true to show observations
      trace1 = scatter3d(;x=[observations[i].y[1] for i in 1:T],
                          y=[observations[i].y[2] for i in 1:T],
                          z=[observations[i].y[3] for i in 1:T],
                          mode="lines",
                          marker=attr(color="#1f77b4", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                          line=attr(color="#1f77b4", width=1))
    else
      #
      trace1 = scatter3d(;x=[truth[i].particles[1,1] for i in 1:T],
                          y=[truth[i].particles[2,1] for i in 1:T],
                          z=[truth[i].particles[3,1] for i in 1:T],
                          mode="lines",
                          marker=attr(color="#1f77b4", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                          line=attr(color="#1f77b4", width=1))
    end
    #
    trace2 = scatter3d(;x=[mean(faps[i].p.particles[1,:],faps[i].p.weights[1]) for i in 1:length(faps)],
                        y=[mean(faps[i].p.particles[2,:],faps[i].p.weights[1]) for i in 1:length(faps)],
                        z=[mean(faps[i].p.particles[3,:],faps[i].p.weights[1]) for i in 1:length(faps)],
                        mode="lines",
                          marker=attr(color="#9467bd", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                        line=attr(color="rgb(44, 160, 44)", width=1))
    #
    trace3 = scatter3d(;x=[mean(fips[i].tip.particles[1,:]) for i in 1:length(fips)],
                        y=[mean(fips[i].tip.particles[2,:]) for i in 1:length(fips)],
                        z=[mean(fips[i].tip.particles[3,:]) for i in 1:length(fips)],
                        mode="lines",
                        marker=attr(color="#bcbd22", size=12, symbol="circle",
                                    line=attr(color="rgb(0,0,0)", width=0)),
                        line=attr(color="#bcbd22", width=1))
    #
    layout = Layout(autosize=false, width=500, height=500,
        margin=attr(l=0, r=0, b=0, t=65))
    plot([trace1, trace2,trace3], layout)
end
evograph()


function fuzzGraph()
    #
    if true #true to show observations
      trace1 = scatter3d(;x=[observations[2].y[1] for i in 1:T],
                          y=[observations[2].y[2] for i in 1:T],
                          z=[observations[2].y[3] for i in 1:T],
                          mode="lines",
                          marker=attr(color="#1f77b4", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                          line=attr(color="#1f77b4", width=1))
    else
      #
      trace1 = scatter3d(;x=[truth[2].particles[1,1] for i in 1:T],
                          y=[truth[2].particles[2,1] for i in 1:T],
                          z=[truth[2].particles[3,1] for i in 1:T],
                          mode="lines",
                          marker=attr(color="#1f77b4", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                          line=attr(color="#1f77b4", width=1))
    end
    #
    trace2 = scatter3d(;x=fips[1].means[1,:],
                        y=fips[1].means[2,:],
                        z=fips[1].means[3,:],
                        mode="markers",
                          marker=attr(color="#9467bd", size=12, symbol="circle",
                                      line=attr(color="rgb(0,0,0)", width=0)),
                        line=attr(color="rgb(44, 160, 44)", width=1))
    #
    trace3 = scatter3d(;x=[mean(fips[1].localDists[1,i])[3] for i in 1:size(fips[1].localDists)[2]],
                        y=[mean(fips[1].localDists[2,i])[3] for i in 1:size(fips[1].localDists)[2]],
                        z=[mean(fips[1].localDists[3,i])[3] for i in 1:size(fips[1].localDists)[2]],
                        mode="markers",
                        marker=attr(color="#bcbd22", size=12, symbol="circle",
                                    line=attr(color="rgb(0,0,0)", width=0)),
                        line=attr(color="#bcbd22", width=1))
    #
    layout = Layout(autosize=false, width=500, height=500,
        margin=attr(l=20, r=0, b=0, t=65))
    plot([trace1, trace2], layout)
end
#fuzzGraph()

length(fips)
length(observations)

i = 1
fp = fips[i]
y = fp.obs
mean(fp.lps) # mean of all quasi-f's
mean([mean(fp.lps[:,i,i]) for i in 1:10]) #mean of same-particle quasi-f's
mean([max(fp.lps[1,i,(i+1):(i+69)]...) for i in 1:10]) #mean of best different-particle quasi-f's

fp.obs.y[1:5]





mean(fp.means,dims=2)[1:5]





mean(fp.base,dims=2)[1:5]





mean(fp.tip.particles,dims=2)[1:5]





[mean(fp.base,fp.ws[j],2)[j] for j in 1:5]





size(fp.tip.particles)
size(fp.ws[1])

observations[i].y[1:5]
print(fp.base[1:2,1],"\n",fp.means[1:2,1])
fips[1].tip.filter.z.h

size(fp.base)
dimover = 2
matrix = fp.base
dimover
fp.ws
bkf.reweight!(fp,fp.obs)
myweights = fp.ws[l]
#function mmean(matrix,myweights,dimover)
outy = [x for x in size(matrix)]
outy[dimover] = 1
outy
size(repeat(myweights,outer=outy))

#sum(matrix .* repeat(myweights,outer=outy),dimover) / sum(myweights)
#end



using StatsBase
d = size(fp.base,1)
fp.ws = Vector{ProbabilityWeights}(undef, d)
diffs = fp.base - fp.tip.filter.z.h * repeat(y.y,outer=[1,fp.tip.n])  # fp.tip.f.z.h should probably be eye ?

fp.base[5,10]
y.y[5]
diffs[5,10]



vars = diag(fp.tip.filter.z.r) #assumes fp.tip.f.z.h is eye and ...r is diagonal
#for l = 1:d

l = 1
if false
    forwardProb = zeros(d)
    for ph = 1:fp.tip.n
        forwardProb += exp.(fp.lps[l,:,ph])
    end
else
    forwardProb = 1. #ones(d)
end
wvec = exp.(-diffs[l,:].^2 ./ vars[l] / 2 ) ./ forwardProb # ./ if forwardProb is vector
#fp.ws[l] = ProbabilityWeights(wvec)
#end
fp.obs=y
results = Array{Float64}(undef,5,3)
results[:,1] = fp.base[l,1:5]
results[:,2] = diffs[l,1:5]
results[:,3] = wvec[1:5]
results
mean(fp.base,ProbabilityWeights(wvec),2)[l]


[mean(fp.base),
  mean(fp.base,ProbabilityWeights(wvec),2)[l],
  y.y[l]]

size(mean(fp.tip.particles,ProbabilityWeights(wvec),2))

sum(fp.ws[2])

ft
y
fp1 = bkf.predictUpdate(ft,y,1)
fp10 = bkf.predictUpdate(ft,y,10)
fp100 = bkf.predictUpdate(ft,y,100)
fp = fp1

[mean(faps[i].p.weights[1]) for i in 1:length(faps)]

bkf.obsNoiseDistribution(faps[i].p.filter)

size(faps)

i= 1
fips[1].obs.y[1:4]




faps[1].y.y[1:4]


i = 2

faps[i].y.y[1:4,1]




mean(faps[i].p.particles[1:4,:],faps[i].p.weights[1],2)




mean(faps[i].p.particles[1:4,:],dims=2)




mean(faps[i].p.weights[1])
using StatsBase
using Statistics
using Distributions
pset = faps[i].p
y = faps[i].y
i = 1
toI = i*pset.hoodSize
fromI = toI - pset.hoodSize + 1
ProbabilityWeights(
          pdf(bkf.obsNoiseDistribution(pset.filter, fromI, toI),
              pset.particles[fromI:toI,:] - pset.filter.z.h[fromI:toI,fromI:toI] * repeat(y.y[fromI:toI],outer=[1,size(pset.particles,2)])
          ))

pset.particles[fromI:toI,1] - pset.filter.z.h[fromI:toI,fromI:toI] * repeat(y.y[fromI:toI],outer=[1,1])
pset.filter.z.h[fromI:toI,fromI:toI]* repeat(y.y[fromI:toI],outer=[1,1])
pset.particles[fromI:toI,1]

ProbabilityWeights(
          pdf(bkf.obsNoiseDistribution(pset.filter, fromI, toI),
              pset.particles[fromI:toI,1] - pset.filter.z.h[fromI:toI,fromI:toI] * repeat(y.y[fromI:toI],outer=[1,1])
          ))

(pset.particles[fromI:toI,:] - pset.filter.z.h[fromI:toI,fromI:toI] *
    repeat(y.y[fromI:toI],outer=[1,size(pset.particles,2)]))[1:4,1:4]
size(pdf(bkf.obsNoiseDistribution(pset.filter, fromI, toI),
    pset.particles[fromI:toI,:] - pset.filter.z.h[fromI:toI,fromI:toI] * repeat(y.y[fromI:toI],outer=[1,size(pset.particles,2)])))

fappy = bkf.FrankenStep(faps[1], faps[2].y)
pset = fappy.p

                          mean(faps[1].p.particles[1:4,:],dims=2)
median(fappy.p.particles[1:4,:],dims=2)
4



smallState = fill(8.,9,2)
smallState[9,1] = 8.1
smallState[10,1] = 8.1
newState = bkf.newCenters(fappy.p.filter,smallState)

smallState = faps[1].p.particles[:,1:4]


mean(newState,dims=2)
