#.....csv
MEquiv = 25
easy = false
clones = 1 #apparent dimensions = d*clones
highcor, corgap, othervar, mainvar, vargap, μ = (.8,2,2.,.5,3,.25)
nIterVec = [200,350,700]
basefname = "convergetruth.csv"
if easy
  difficulty = "easy_"
  timeSuperStep = 0.05
else
  difficulty = "hard_"
  timeSuperStep = 0.4
end
fname = difficulty*basefname
outprefix = "outcome_converged2_"
outcomefile = outprefix*difficulty*string(MEquiv)*"_"*ENV["USER"]*".csv"
full_d = 20 #dimensions
s = 2 #steps
d = div(full_d, clones)
doAlgos = true

function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

kalmandir = join(split(Base.source_path(),'/')[1:end-2],'/')
ppath(kalmandir * "/src")
if false #kalmandir == "/Users/chema/Dropbox/Kalman.jl"
  ppath("/Users/chema/Dropbox/")
  ppath("/Users/chema/mydev/Gadfly.jl/src")
end



using Revise
using Distributions, DataStructures
using bkf
using DelimitedFiles
using Random
using LinearAlgebra


mymodel = bkf.createLorenzModel(d, timeSuperStep)
bkf.makeFunkyInitialDist!(mymodel, highcor, corgap, othervar, mainvar, vargap, μ)

mydict = OrderedDict()

kfa = bkf.KfAlgo()
bkf.putParams!(kfa,mydict)
bkf.init(kfa, mymodel)

pfa = bkf.PfAlgo(MEquiv)
bkf.putParams!(pfa, mydict)
bkf.init(pfa, mymodel)

ba = bkf.BlockAlgo(MEquiv,4)
bkf.putParams!(ba, mydict)
bkf.init(ba, mymodel)


fa30u = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,
                    30, #histPerLoc
                    100, #nIter
                    1., #useForward
                    4, #overlap
                    bkf.FuzzFinkelParticles,
                    .125, #rejuv
                    )
#
fa10u = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,
                    10, #histPerLoc
                    100, #nIter
                    1., #useForward
                    4, #overlap
                    bkf.FuzzFinkelParticles,
                    .125, #rejuv
                    )
#
fa30l = bkf.FinkelAlgo(MEquiv,1,bkf.SampleLog,bkf.MhSampled,
                    30, #histPerLoc
                    100, #nIter
                    1., #useForward
                    4, #overlap
                    bkf.FuzzFinkelParticles,
                    .125, #rejuv
                    )
#
fa10l = bkf.FinkelAlgo(MEquiv,1,bkf.SampleLog,bkf.MhSampled,
                    10, #histPerLoc
                    100, #nIter
                    1., #useForward
                    4, #overlap
                    bkf.FuzzFinkelParticles,
                    .125, #rejuv
                    )
#
bkf.putParams!(fa30u, mydict)
#bkf.init(fa, mymodel)

obs = 0
try
    global obs = bkf.loadObservations(fname)
    print(obs[3][3].x.x[3])
    print("QQQQQQQQQQQQQQQQQQQ")
catch
    #@assert "Don't recreate; too late." == 0
    global obs = bkf.createObservations(mymodel, s)
    bkf.saveObservations(obs, fname, false, clones)
    if clones>1
      bkf.saveObservations(obs, "uncloned_"*fname, false, 1)
    end
end

algos = vcat([ba],bkf.finkelAlgos(MEquiv))

if doAlgos
  bkf.testConvergence(mymodel, [fa30u,fa10u,fa30l,fa10u], nIterVec, 360, outcomefile)
end
