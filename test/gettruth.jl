
MEquiv = 40
easy = false
clones = 1 #apparent dimensions = d*clones
basefname = "truth.csv"
if easy
  difficulty = "easy_"
  timeSuperStep = 0.05
else
  difficulty = "hard_"
  timeSuperStep = 0.4
end
fname = "uncloned_"*difficulty*basefname
outcomefile = "truth_"*difficulty*string(MEquiv)*".csv"
full_d = 5 #dimensions
s = 60 #steps
d = div(full_d, clones)
if clones>1
  fname = "repeating_" * fname
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


mod = bkf.createLorenzModel(d, timeSuperStep)

mydict = OrderedDict()

kfa = bkf.KfAlgo()
bkf.putParams!(kfa,mydict)
bkf.init(kfa, mod)

pf = bkf.PfAlgo(MEquiv)
bkf.putParams!(pf, mydict)
bkf.init(pf, mod)

ba = bkf.BlockAlgo(MEquiv,4)
bkf.putParams!(ba, mydict)
bkf.init(ba, mod)

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
                    MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS, #overlap
                    bkf.FuzzFinkelParticles,
                    1/(MEquiv^(1-1/bkf.DEFAULT_PRODUCT_RADIUS)/bkf.DEFAULT_PRODUCT_RADIUS) /2, #rejuv
                    )
#
bkf.putParams!(faLog, mydict)
bkf.putParams!(fa1050, mydict)
#bkf.init(fa, mod)



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
  bkf.runAlgos(mod, obs, [pf], 1, outcomefile)
end
