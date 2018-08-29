using DataFrames, CSV, DataStructures

##############filtering algorithms

abstract type Algo end #filtering algorithms

function predictupdate(state, obs, algo::Algo) #override if you need to pass in extra params - ie, nIter
    predictupdate(state, obs)
end

function resample(state, algo::Algo)
    s = state
    g = "(no)\n"
    try
        s = resample(state)
        g = "(yes)\n"
        print("Resampled\n")
    catch y
        print("Couldn't resample: ",typeof(y),"\n")
    end
    print(g)
    s
end

type KfAlgo <: Algo
end

function putParams!(algo::KfAlgo,row::OrderedDict)
    row[:model] = "Kalman"
    row[:mod] = "kal"
    row
end

function init(algo::KfAlgo, model::KalmanFilter)
    model
end

type PfAlgo <: Algo
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

type BlockAlgo <: Algo
    MEquiv::Int64
    r::Int64
end

function getM(algo::BlockAlgo)
    div(algo.MEquiv^2, 5)
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
    FrankenStep(model, getM(algo), algo.r)
end

type FinkelAlgo <: Algo
    MEquiv::Int64
    r::Int64
    sampType::Type #for now, fixed subparams
    mhType::Type
    histPerLoc::Int64
    #rsub::Int64 #possibly meaningless
    nIter::Int64
    useForward::Float64
end

function predictupdate(state, obs, algo::FinkelAlgo) #override if you need to pass in extra params - ie, nIter
    predictupdate(state, obs, algo.nIter)
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
    row
end

function init(algo::FinkelAlgo, model::KalmanFilter)
    FinkelToe(model, getM(algo), FinkelParams(algo.sampType(),
                                        algo.mhType(algo.r,algo.histPerLoc),
                                        algo.useForward))
end

#TODO: algos for kalman "idealized" versions of block and finkel

function createModel(d)
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

    x0 = State(zeros(d),initialvar*eye(d))

    a = Tridiagonal(bleedl * ones(d-1),bleedm * ones(d),bleedr * ones(d-1))
    a = a * temper / (bleedl+bleedm+bleedr) #progression matrix
    #a = a * temper / det(a)^(1/d) #progression matrix - divergent

    jitvec = fill(lojitter,d)
    jitvec[1:jitgap:d] = hijitter
    g = Diagonal(jitvec)


    q = eye(d)

    f = LinearModel(a,g,q)


    noisevec = fill(basenoise,d)
    noisevec[1:lowgap:d] = lownoise
    r = Diagonal(noisevec)
    #var35 = meanvarlocs(zeros(d), inv(h), 3:5)[2]
    z = LinearObservationModel(Array(r))

    kf0 = BasicKalmanFilter(x0,f,z)

    kf0
end

function createObservations(model, steps)
    pf1 = ParticleSet(model,1)

    truth = Vector{ParticleSet}(0)#length(t))
    observations = Vector{Observation}(0)#length(t))
    kfs = Vector{BasicKalmanFilter}(0)#length(t))
    kf = model

    push!(kfs, kf)

    push!(truth, pf1)
    push!(observations, Observation(pf1,1))

    for i in 2:(steps+1)

        push!(truth, ap(truth[i-1]))
        push!(observations, Observation(truth[i],1))
        kf2 = bkf.predict(kf)
        kf = update(kf2,observations[i])
        push!(kfs, kf)
    end

    (truth, observations, kfs)
end

function trsp(v)
    reshape(v,(1,:))
end

function saveObservations(obs, fileName, overwrite = false)
    if !overwrite
        assert(!isfile(fileName))
    elseif isfile(fileName)
        rm(fileName)
    end
    (truth, observations, kfs) = obs
    for i=2:length(truth)
        open( fileName,  "a") do outfile
            writecsv( outfile, trsp([truth[i].particles[:,1]
                                     observations[i].y
                                    ]))
        end
    end
end

function loadObservations(fileName)
    data = readcsv(fileName)
    d = div(size(data)[2],2)
    model = createModel(d)

    pf1 = ParticleSet(model,1)

    truth = Vector{ParticleSet}(0)#length(t))
    observations = Vector{Observation}(0)#length(t))
    kfs = Vector{BasicKalmanFilter}(0)#length(t))
    kf = model

    push!(kfs, kf)

    push!(truth, pf1)
    push!(observations, Observation(pf1,1)) #junk data - not used

    for i in 1:(size(data)[1])

        push!(truth, ParticleSet(model,
                                1,
                                data[i:i,1:d],
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
        print("Appending to existing save file! \n")
        #check file validity
    else
        print("Writing column headers\n")
        names = [string(param) for param in keys(blankParams)]
        push!(names,"dimension")
        push!(names,"rep")
        push!(names,"startTime")
        push!(names,"step")
        push!(names,"checksum")
        push!(names,"putime")
        push!(names,"samptime")
        push!(names,"klerror")
        names = vcat(names,["kl","covdiv","meandiv","entropydiv"])
        names = vcat(names,[("b"*lpad(i,2,0)) for i in 1:d])
        names = vcat(names,[("v"*lpad(i,2,0)^2) for i in 1:d])
        names = vcat(names,[("v"*lpad(i,2,0)*lpad(i+1,2,0)) for i in 1:(d-1)])
        names = vcat(names,[("v"*lpad(i,2,0)*lpad(i+2,2,0)) for i in 1:(d-2)])

        open( saveFileName,  "a") do outfile
            writecsv( outfile, trsp(names))
        end
    end




    for rep in 1:reps
        for algo = algos
            paramDict = deepcopy(blankParams)
            putParams!(algo,paramDict)
            basedatavec = [get(paramDict,k,"") for k in keys(paramDict)]
            print("\n\nAlgo: ",basedatavec,"\n")
            push!(basedatavec,d) #dimension
            push!(basedatavec,rep) #rep
            push!(basedatavec,runstarttime) #time

            state = init(algo, model)
            for i = 2:length(truth)
                print("Step ",i-1,":     ",basedatavec,"\n")
                putime = (@timed state = predictupdate(state, observations[i], algo))[2]
                #measure something here? maybe TODO later
                samptime = (@timed state = resample(state, algo))[2]
                (μ2,Σ2) = musig(state)
                (μ1,Σ1) = musig(kfs[i])
                μ = μ2 - μ1

                datavec = copy(basedatavec)
                push!(datavec,i-1) #step
                push!(datavec,kfs[i].x.x[d]) #checksum
                push!(datavec,putime) #putime
                push!(datavec,samptime) #samptime
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
                                [Σ2[i,i] for i in 1:d],
                                [Σ2[i,i+1] for i in 1:(d-1)],
                                [Σ2[i,i+2] for i in 1:(d-2)])

                open( saveFileName,  "a") do outfile
                    writecsv( outfile, trsp(datavec))
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
                                                    useForward
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
