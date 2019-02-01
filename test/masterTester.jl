MEquiv = 4500
easy = true
hardMeansHard = false
useRepeats = false
clones = 1 #apparent dimensions = d*clones
doTestConvergence = false
basefname = "truth.csv"
outprefix = "outcome_"
full_d = 5 #dimensions
s = 60 #steps
doAlgos = true
neighborhood_r = 4
overlap = exp(.5)*2
rejuv = 1/overlap/4
histPerLoc = 15
nIter = 150
nIterVec = [1,30,100,250,600,1300]

#setup for imports
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


#imports
using Revise
using Distributions, DataStructures
using bkf
using DelimitedFiles
using Random
using LinearAlgebra

#use global variables to set up file names
if easy
  difficulty = "easy_"
  timeSuperStep = 0.05
else
  if hardMeansHard
    difficulty = "nightmare_"
    timeSuperStep = 0.6
  else
    difficulty = "hard_"
    timeSuperStep = 0.4
  end
end
convPart = doTestConvergence ? "converge_" : ""
fname = difficulty*basefname
outcomefile = outprefix*string(bkf.REPEATABLE_VERSION)*convPart*difficulty*string(MEquiv)*"_"*ENV["USER"]

if useRepeats #clones>1
  d = div(full_d, clones)
  fname = "repeating_" * fname
else
  d = full_d
  outcomefile = outcomefile*"_nonrep"
end
outcomefile = outcomefile*".csv"
if MEquiv >= 3125
  fname = "uncloned_repeating_easy_truth.csv"
end


#setup model and algo variables
mymodel = bkf.createLorenzModel(d, timeSuperStep)

mydict = OrderedDict()

kfa = bkf.KfAlgo()
bkf.putParams!(kfa,mydict)
bkf.init(kfa, mymodel)

pfa = bkf.PfAlgo(MEquiv)
bkf.putParams!(pfa, mydict)
bkf.init(pfa, mymodel)

ba = bkf.BlockAlgo(MEquiv,neighborhood_r)
bkf.putParams!(ba, mydict)
bkf.init(ba, mymodel)

fa = bkf.FinkelAlgo(MEquiv,neighborhood_r,
                    bkf.SampleLog,bkf.MhSampled,
                    histPerLoc, #histPerLoc
                    nIter, #nIter
                    1., #useForward
                    overlap, #overlap
                    bkf.FuzzFinkelParticles,
                    rejuv, #rejuv
                    )
#
faUni = bkf.FinkelAlgo(MEquiv,neighborhood_r,
                    bkf.SampleUniform,bkf.MhSampled,
                    histPerLoc, #histPerLoc
                    nIter, #nIter
                    1., #useForward
                    overlap, #overlap
                    bkf.FuzzFinkelParticles,
                    rejuv, #rejuv
                    )
#
bkf.putParams!(fa, mydict)
#bkf.init(fa, mymodel)

obs = 0
try
    global obs = bkf.loadObservations(fname)
    print(obs[3][3].x.x[3])
    print("QQQQQQQQQQQQQQQQQQQ")
catch
    @assert "Don't recreate; too late." == 0
    global obs = bkf.createObservations(mymodel, s)
    bkf.saveObservations(obs, fname, false, clones)
    if clones>1
      bkf.saveObservations(obs, "uncloned_"*fname, false, 1)
    end
end

algos = vcat([ba],bkf.finkelAlgos(MEquiv))

#run!
if doAlgos
  if doTestConvergence
    bkf.testConvergence(mymodel, [fa,faUni], nIterVec, 360, outcomefile)
  else
    bkf.runAlgos(mymodel, obs, [pfa], 360, outcomefile)
  end
end
