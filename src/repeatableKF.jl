using DataFrames, CSV, DataStructures, Dates
using NaNMath

REPEATABLE_VERSION = 1.3

##############filtering algorithms

myWritecsv(io, a; opts...) = writedlm(io, a, ','; opts...)
myReadcsv(io; opts...) = readdlm(io, ','; opts...)

abstract type Algo end #filtering algorithms

function predictUpdate(state, obs, algo::Algo) #override if you need to pass in extra params - ie, nIter
    predictUpdate(state, obs)
end


mutable struct KfAlgo <: Algo
end

function putParams!(algo::KfAlgo,row::OrderedDict)
    row[:model] = "Kalman"
    row[:mod] = "kal"
    row
end

function init(algo::KfAlgo, model::KalmanFilter)
    model
end

mutable struct PfAlgo <: Algo
    MEquiv::Int64
end

function getM(algo::PfAlgo)
    algo.MEquiv ^2
end

function putParams!(algo::PfAlgo,row::OrderedDict)
    row[:model] = "Bootstrap PF"
    row[:mod] = "boot"
    row[:M] = getM(algo)
    row[:MEquiv] = algo.MEquiv
    row
end

function init(algo::PfAlgo, model::KalmanFilter)
    ParticleStep(model,getM(algo))
end

mutable struct BlockAlgo <: Algo
    MEquiv::Int64
    r::Int64
end

function getM(algo::BlockAlgo)
    div(algo.MEquiv^2, 2)
end

function putParams!(algo::BlockAlgo,row::OrderedDict)
    row[:model] = "Block PF"
    row[:mod] = "block"
    row[:M] = getM(algo)
    row[:MEquiv] = algo.MEquiv
    row[:r] = algo.r
    row
end

function init(algo::BlockAlgo, model::KalmanFilter)
    toFrankenSet(model, getM(algo), algo.r) #kinda a hack
end

mutable struct FinkelAlgo <: Algo
    MEquiv::Int64
    r::Int64
    sampType::Type #for now, fixed subparams
    mhType::Type
    histPerLoc::Int64
    #rsub::Int64 #possibly meaningless
    nIter::Int64
    useForward::Float64
    overlap::Float64
    algo::Type
    rejuv::Float64
end

function predictUpdate(state, obs, algo::FinkelAlgo) #override if you need to pass in extra params - ie, nIter
    predictUpdate(state, obs, algo.nIter)
end

function getM(algo::FinkelAlgo)
    algo.MEquiv
end

function putParams!(algo::FinkelAlgo,row::OrderedDict)
    row[:model] = "Finkelstein PF"
    row[:mod] = "finkel"
    row[:M] = getM(algo)
    row[:MEquiv] = algo.MEquiv
    row[:sampType] = string(algo.sampType)
    row[:mhType] = string(algo.mhType)
    row[:hpl] = algo.histPerLoc
    row[:r] = algo.r
    #row[:rsub] = algo.rsub
    row[:nIter] = algo.nIter
    row[:useForward] = algo.useForward
    row[:overlap] = algo.overlap
    row[:isFuzzed] = string(algo.algo)
    row[:rejuv] = algo.rejuv
end

function init(algo::FinkelAlgo, model::KalmanFilter)
    FinkelToe(model, getM(algo), FinkelParams(algo.sampType(),
                                        algo.mhType(algo.r,algo.histPerLoc),
                                        algo.useForward,
                                        algo.overlap,
                                        algo.algo,
                                        algo.rejuv
                                        ))
end

#TODO: algos for kalman "idealized" versions of block and finkel

function createLinearModel(d)
    bleedl = .7#1
    bleedm = .8
    bleedr = .1#.25
    hijitter = 1.
    lojitter = .25
    jitgap = 2
    temper = .8 #1/sqrt(2)
    basenoise = 1.
    lownoise = .16
    lowgap = 5
    noisebleed = 0.
    mciter = 6
    initialvar = 1/(1-temper)

    x0 = State(zeros(d),initialvar*Matrix(1.0I,d,d))

    a = Tridiagonal(bleedl * ones(d-1),bleedm * ones(d),bleedr * ones(d-1))
    a = a * temper / (bleedl+bleedm+bleedr) #progression matrix
    #a = a * temper / det(a)^(1/d) #progression matrix - divergent

    jitvec = fill(lojitter,d)
    jitvec[1:jitgap:d] = hijitter
    g = Diagonal(jitvec)


    q = Matrix(1.0I,d,d)

    f = LinearModel(a,g,q)


    noisevec = fill(basenoise,d)
    noisevec[1:lowgap:d] = lownoise
    r = Diagonal(noisevec)
    #var35 = meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = LinearObservationModel(Array(r))

    kf0 = BasicKalmanFilter(x0,f,z)

    kf0
end


function createLorenzModel(d, timeSuperStep = 0.05)

    forcingF = 8.
    timeStep = 0.002 #timeSuperStep/numSteps
    numSteps = div(timeSuperStep,timeStep)
    @assert numSteps > 20
    processNoiseVar = 1e-100
    measurementNoiseVar = [1.]
    mnvvec = repeat(measurementNoiseVar,40)[1:d]
    initialvar = 0.09



    x0 = bkf.State(ones(d)*forcingF,initialvar*Matrix(1.0I,d,d))

    g = Matrix(1.0I,d,d) * processNoiseVar

    q = Matrix(1.0I,d,d)

    f = bkf.LorenzModel(forcingF,d,timeStep,numSteps,g,q)#,overlapFactor)

    r = diagm(0=>mnvvec)

    #var35 = bkf.meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = bkf.LinearObservationModel(Array(r))
    var35 = bkf.meanvarlocs(zeros(d), inv(r), 3:5)[2]

    kf0 = bkf.BasicLorenzFilter(x0,f,z)
    kf0
end

function createObservations(model, steps)
    pf1 = ParticleSet(model,1)

    truth = Vector{ParticleSet}(undef,0)#length(t))
    observations = Vector{Observation}(undef,0)#length(t))
    kfs = Vector{KalmanFilter}(undef,0)#length(t))
    difs = Vector{Vector{Float64}}(undef,0)
    kf = model

    push!(kfs, kf)

    push!(truth, pf1)
    push!(observations, Observation(pf1,1))
    push!(difs, [observations[1].y[l]-truth[1].particles[l,1] for l in 1:5])


    for i in 2:(steps+1)

        push!(truth, ap(truth[i-1]))
        push!(observations, Observation(truth[i],1))
        push!(difs, [observations[i].y[l]-truth[i].particles[l,1] for l in 1:5])
        #debug(difs[i])
        kf2 = bkf.predict(kf)
        kf = update(kf2,observations[i])
        push!(kfs, kf)
    end

    (truth, observations, kfs)
end

if false #testing code; pseudo-pastebin
  model = bkf.createLorenzModel(20)
  pf1 = bkf.ParticleSet(model,1)
  (mytruths, myobss, mykfs,mydifs) = bkf.createObservations(model,100)
  [(mean(myobss[j].y[l]-mytruths[j].particles[l,1] for j in 2:100),
    var(myobss[j].y[l]-mytruths[j].particles[l,1] for j in 2:100),
    cor([myobss[j].y[l] for j in 2:100],[mytruths[j].particles[l,1] for j in 2:100]))
          for l in 1:20]
  l = 1
  [(myobss[j].y[l]-mytruths[j].particles[l,1],mydifs[j][l]) for j in 1:5]
end

function trsp(v)
    reshape(v,(1,:))
end

function saveObservations(obs, fileName, overwrite = false, clones = 1)
    if !overwrite
        @assert(!isfile(fileName))
    elseif isfile(fileName)
        rm(fileName)
    end
    (truth, observations, kfs) = obs
    for i=1:length(truth)
        open( fileName,  "a") do outfile
            myWritecsv( outfile, trsp([repeat(truth[i].particles[:,1],clones)
                                     repeat(observations[i].y,clones)
                                    ]))
        end
    end
end

function loadObservations(fileName)
    data = myReadcsv(fileName)
    d = div(size(data)[2],2)
    model = createLorenzModel(d)

    pf1 = ParticleSet(model,1)

    truth = Vector{ParticleSet}(undef,0)#length(t))
    observations = Vector{Observation}(undef,0)#length(t))
    kfs = Vector{KalmanFilter}(undef,0)#length(t))
    kf = model

    push!(kfs, kf)

    push!(truth, pf1)
    push!(observations, Observation(pf1,1)) #junk data - not used

    for i in 1:(size(data)[1])

        push!(truth, ParticleSet(model,
                                1,
                                Matrix((data[i:i,1:d])'),
                                ProbabilityWeights(ones(d))))

        push!(observations,
            Observation(
                vec(
                    data[i,((d+1):(2*d))]
                    )))
        kf2 = bkf.predict(kf)
        kf = update(kf2,observations[i+1])
        push!(kfs, kf)
    end

    (truth, observations, kfs)
end

function runAlgos(model, obs, algos, reps, saveFileName)
    #debug(obs)
    (truth, observations, kfs) = obs
    runstarttime = Dates.now()
    d = length(model.x.x)
    blankParams = OrderedDict()
    for algo = algos
        putParams!(algo,blankParams)
    end
    i = 1
    for k in keys(blankParams)
        blankParams[k] = ""
    end

    if isfile(saveFileName)
        print("Appending to existing save file! \n"); print("\n","""print("Appending to existing save file! \n")""")
        #check file validity
    else
        print("Writing column headers\n")
        names = [string(param) for param in keys(blankParams)]
        push!(names,"version")
        push!(names,"dimension")
        push!(names,"rep")
        push!(names,"startTime")
        push!(names,"step")
        push!(names,"checksum")
        push!(names,"putime")
        push!(names,"samptime")
        push!(names,"klerror")
        names = vcat(names,["kl","covdiv","meandiv","entropydiv"])
        names = vcat(names,[("d"*lpad(i,2,"0")) for i in 1:d])
        names = vcat(names,[("v"*lpad(i,2,"0")^2) for i in 1:d])
        names = vcat(names,[("v"*lpad(i,2,"0")*lpad(i+1,2,"0")) for i in 1:(d-1)])
        names = vcat(names,[("v"*lpad(i,2,"0")*lpad(i+2,2,"0")) for i in 1:(d-2)])

        names = vcat(names,[("t"*lpad(i,2,"0")) for i in 1:d])
        names = vcat(names,[("o"*lpad(i,2,"0")) for i in 1:d])
        open( saveFileName,  "a") do outfile
            myWritecsv( outfile, trsp(names))
        end
    end




    for rep in 1:reps
        debug("Rep:",rep)
        debug()
        debug()
        debug()
        for algo = algos
            paramDict = deepcopy(blankParams)
            putParams!(algo,paramDict)
            basedatavec = [get(paramDict,k,"") for k in keys(paramDict)]
            print("\n\nAlgo: ",basedatavec,"\n"); print("\n","""print("\n\nAlgo: ",basedatavec,"\n")""")
            push!(basedatavec,REPEATABLE_VERSION) #version
            push!(basedatavec,d) #dimension
            push!(basedatavec,rep) #rep
            push!(basedatavec,runstarttime) #time

            state = init(algo, model)
            for i = 2:length(truth)
                print("Step ",i-1,":     ",basedatavec,"\n"); print("\n","""print("Step ",i-1,":     ",basedatavec,"\n")""")
                putime = (@timed state = predictUpdate(state, observations[i], algo))[2]
                #measure something here? maybe TODO later
                (μ2,Σ2) = musig(state)
                (μ1,Σ1) = musig(kfs[i])
                debug("sizes",i,size(μ2),size(truth[i].particles))
                μ = μ2 - truth[i].particles[:,1]
                print("Meansqs: resampled:",mean(μ.^2),"\n"); print("\n","""print("Meansqs: resampled:",mean(μ.^2),"\n")""")

                datavec = copy(basedatavec)
                push!(datavec,i-1) #step
                push!(datavec,kfs[i].x.x[d]) #checksum
                push!(datavec,putime) #putime
                push!(datavec,0) #samptime
                try
                    kl = kl2(μ1,Σ1,μ2,Σ2)
                    push!(datavec,"") #no error
                    append!(datavec,kl)
                catch y
                    push!(datavec,string(typeof(y)))
                    append!(datavec,["","","",""])
                end
                datavec = vcat(datavec,
                                μ,
                                [Σ2[l,l] for l in 1:d],
                                [Σ2[l,l+1] for l in 1:(d-1)],
                                [Σ2[l,l+2] for l in 1:(d-2)],
                                [truth[i].particles[l,1] for l in 1:d],
                                [observations[i].y[l] for l in 1:d]
                              )

                open( saveFileName,  "a") do outfile
                    myWritecsv( outfile, trsp(datavec))
                end
            end
        end
    end
end

function finkelAlgos(Ms, sampTypes=[SampleUniform], mhTypes=[MhSampled], histPerLocs=[30], nIters=[160, 0], useForwards=[1.])
    algos = Algo[]
    for samp in sampTypes
        for mh in mhTypes
            for histPerLoc in histPerLocs
                for useForward in useForwards
                    for nIter in nIters
                        for M in Ms
                            push!(algos, FinkelAlgo(M,
                                                    1, #default r. TODO, fix
                                                    samp,
                                                    mh,
                                                    histPerLoc,
                                                    nIter,
                                                    useForward,
                                                    2., FuzzFinkelParticles, .25
                                                    )
                                )
                        end
                    end
                end
            end
        end
    end
    algos
end
