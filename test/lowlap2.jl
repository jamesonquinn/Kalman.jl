#outcome_lowlap_2_hard_250_nonrep.csv
MEquiv = 250
easy = false
useRepeats = false
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
outprefix = "outcome_lowlap_2_"
outcomefile = outprefix*difficulty*string(MEquiv)*".csv"
full_d = 40 #dimensions
s = 60 #steps
d = div(full_d, clones)
if useRepeats #clones>1
  fname = "repeating_" * fname
else
  outcomefile = outprefix*difficulty*string(MEquiv)*"_nonrep.csv"
end
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

faLog = bkf.FinkelAlgo(MEquiv,1,bkf.SampleLog,bkf.MhSampled,
                    30, #histPerLoc
                    60, #nIter
                    1., #useForward
                    MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS, #overlap
                    bkf.FuzzFinkelParticles,
                    1/(MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS) /2, #rejuv
                    )
#
fa1050 = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,
                    10, #histPerLoc
                    50, #nIter
                    1., #useForward
                    4, #overlap
                    bkf.FuzzFinkelParticles,
                    .125, #rejuv
                    )
#
bkf.putParams!(faLog, mydict)
bkf.putParams!(fa1050, mydict)
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
  bkf.runAlgos(mymodel, obs, [fa1050], 360, outcomefile)
end
