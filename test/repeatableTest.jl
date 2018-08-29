function ppath(p)
  if LOAD_PATH[end] != p
    push!(LOAD_PATH,p)
  else
    println(p," already in LOAD_PATH.")

  end
end

ppath("/Users/chema/mydev/Kalman.jl/src")
ppath("/Users/chema/mydev/")
ppath("/Users/chema/mydev/Gadfly.jl/src")

macro load(pkg)
  quote
    if isdefined(Symbol($(string(pkg))))
      #try
      #  reload(ucfirst(string(pkg)))
      #catch
        reload($(string(pkg)))
      #end
    else
      import $pkg
    end
  end
end

using Distributions, DataStructures

@load bkf



d = 30
s = 10
MEquiv = 400

mod = bkf.createModel(d)

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

fa = bkf.FinkelAlgo(MEquiv,1,bkf.SampleUniform,bkf.MhSampled,30,100,1.)
bkf.putParams!(fa, mydict)
bkf.init(fa, mod)


obs = bkf.createObservations(mod, 5)

testfname = "dummyworld.csv"
bkf.saveObservations(obs, testfname, true)

obs2 = bkf.loadObservations(testfname)
obs2[3][6].x.x .- obs[3][6].x.x #should be zeros(mydict)

fname = "myworld.csv"
obs = 0
if false
    obs = bkf.createObservations(mod, s)
    bkf.saveObservations(obs, fname, false)
else
    obs = bkf.loadObservations(fname)
    print(obs[3][3].x.x[3])
end

algos = vcat([ba],bkf.finkelAlgos(MEquiv))
outcomefile = "outcomes_"*string(MEquiv)*".csv"

bkf.runAlgos(mod, obs, algos, 30, outcomefile)
