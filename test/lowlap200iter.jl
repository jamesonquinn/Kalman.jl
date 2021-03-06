#outcome_lowlap_2_hard_250_nonrep.csv
MEquiv = 250
easy = false
useRepeats = true
clones = 1 #apparent dimensions = d*clones
basefname = "truth.csv"
if easy
  difficulty = "easy_"
  timeSuperStep = 0.05
else
  difficulty = "hard_"
  timeSuperStep = 0.4
end
fname = difficulty*basefname
outprefix = "outcome_250iter_"
outcomefile = outprefix*difficulty*string(MEquiv)*"_"*ENV["USER"]*".csv"
full_d = 40 #dimensions
s = 60 #steps
d = div(full_d, clones)
if useRepeats #clones>1
  fname = "repeating_" * fname
else
  outcomefile = outprefix*difficulty*string(MEquiv)*"_nonrep.csv"
end
doAlgos = true
neighborhood_r = 4
overlap = exp(.5)*2
rejuv = 1/overlap/4
histPerLoc = 15
nIters = 250
nIterVec = [1,30,100,250,600,1300]
doTestConvergence = true 

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
                    15, #histPerLoc
                    200, #nIter
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

if doAlgos
  bkf.runAlgos(mymodel, obs, [fa], 360, outcomefile)
end
