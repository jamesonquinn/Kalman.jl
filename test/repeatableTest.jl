
MEquiv = 261
easy = false
clones = 1 #apparent dimensions = d*clones
if easy
  fname = "easy_truth.csv"
  timeSuperStep = 0.05
else
  fname = "hard_truth.csv"
  timeSuperStep = 0.4
end
outcomefile = "myworld3__"*string(MEquiv)*".csv"
full_d = 40 #dimensions
s = 60 #steps
d = div(full_d, clones)
if clones>1
  fname = "repeating_" * fname
end
doAlgos = false

function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

ppath("/Users/chema/Dropbox/Kalman.jl/src")
ppath("/Users/chema/Dropbox/")
ppath("/Users/chema/mydev/Gadfly.jl/src")



using Revise
using Distributions, DataStructures
using bkf
using DelimitedFiles
using Random
using LinearAlgebra


mod = bkf.createLorenzModel(d, timeSuperStep)

mydict = OrderedDict()

kfa = bkf.KfAlgo()
bkf.putParams!(kfa,mydict)
bkf.init(kfa, mod)

pfa = bkf.PfAlgo(MEquiv)
bkf.putParams!(pfa, mydict)
bkf.init(pfa, mod)

ba = bkf.BlockAlgo(MEquiv,5)
bkf.putParams!(ba, mydict)
bkf.init(ba, mod)

fa1 = bkf.FinkelAlgo(MEquiv,1,bkf.SampleLog,bkf.MhSampled,
                    30, #histPerLoc
                    60, #nIter
                    1., #useForward
                    MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS, #overlap
                    bkf.FuzzFinkelParticles,
                    1/(MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS) /2, #rejuv
                    )
#
fa2 = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,
                    15, #histPerLoc
                    100, #nIter
                    1., #useForward
                    MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS, #overlap
                    bkf.FuzzFinkelParticles,
                    1/(MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS) /2, #rejuv
                    )
#
bkf.putParams!(fa1, mydict)
bkf.putParams!(fa2, mydict)
#bkf.init(fa, mod)


obs = bkf.createObservations(mod, 5)

testfname = "dummy2world.csv"
bkf.saveObservations(obs, testfname, true)

obs2 = bkf.loadObservations(testfname)
obs2[3][6].x.x .- obs[3][6].x.x #should be zeros(mydict)

obs = 0
try
    global obs = bkf.loadObservations(fname)
    print(obs[3][3].x.x[3])
    print("QQQQQQQQQQQQQQQQQQQ")
catch
    global obs = bkf.createObservations(mod, s)
    bkf.saveObservations(obs, fname, false, clones)
    if clones>1
      bkf.saveObservations(obs, "uncloned_"*fname, false, 1)
    end
end

algos = vcat([ba],bkf.finkelAlgos(MEquiv))

if doAlgos
  bkf.runAlgos(mod, obs, [fa1, fa2], 360, outcomefile)
end
