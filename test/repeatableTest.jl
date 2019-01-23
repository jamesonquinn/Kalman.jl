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

d = 20
s = 10
MEquiv = 40

mod = bkf.createLorenzModel(d)

mydict = OrderedDict()

kfa = bkf.KfAlgo()
bkf.putParams!(kfa,mydict)
bkf.init(kfa, mod)

pfa = bkf.PfAlgo(MEquiv)
bkf.putParams!(pfa, mydict)
bkf.init(pfa, mod)

ba = bkf.BlockAlgo(MEquiv,3)
bkf.putParams!(ba, mydict)
bkf.init(ba, mod)

fa = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,
                    30, #histPerLoc
                    40, #nIter
                    1., #useForward
                    2., #overlap
                    bkf.FuzzFinkelParticles,
                    .25, #rejuv
                    )
bkf.putParams!(fa, mydict)
bkf.init(fa, mod)


obs = bkf.createObservations(mod, 5)

testfname = "dummy2world.csv"
bkf.saveObservations(obs, testfname, true)

obs2 = bkf.loadObservations(testfname)
obs2[3][6].x.x .- obs[3][6].x.x #should be zeros(mydict)

fname = "myLorenzWorld.csv"
obs = 0
if false
    obs = bkf.createObservations(mod, s)
    bkf.saveObservations(obs, fname, false)
else
    obs = bkf.loadObservations(fname)
    print(obs[3][3].x.x[3]); print("
","""print(obs[3][3].x.x[3])""")
end

algos = vcat([ba],bkf.finkelAlgos(MEquiv))
outcomefile = "outcomes_big_"*string(MEquiv)*".csv"

bkf.runAlgos(mod, obs, algos, 3, outcomefile)
