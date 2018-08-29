using DataFrames, CSV

##############filtering algorithms

abstract type Algo end #filtering algorithms

type KfAlgo <: Algo
end

function putParams!(algo::KfAlgo,row::Dict)
    row[:model] = "Kalman"
    row[:mod] = "kal"
end

function init(algo::KfAlgo, model::KalmanFilter)
    model
end

type PfAlgo <: Algo
    M::Int64
end

function putParams!(algo::PfAlgo,row::Dict)
    row[:model] = "Bootstrap PF"
    row[:mod] = "boot"
    row[:M] = algo.M
    row[:MEquiv] = round(Int64,sqrt(algo.M))
end

function init(algo::PfAlgo, model::KalmanFilter)
    ParticleStep(model,algo.M)
end

type BlockAlgo <: Algo
    M::Int64
    r::Int64
end

function putParams!(algo::BlockAlgo,row::Dict)
    row[:model] = "Block PF"
    row[:mod] = "block"
    row[:M] = algo.M
    row[:MEquiv] = round(Int64,sqrt(algo.M * 5))
    row[:r] = algo.r
end

function init(algo::BlockAlgo, model::KalmanFilter)
    FrankenStep(model, algo.M, algo.r)
end

type FinkelAlgo <: Algo
    M::Int64
    sampType::Type #for now, fixed subparams
    mhType::Type
    histPerLoc::Int64
    r::Int64
    rsub::Int64 #possibly meaningless
    nIter::Int64
    #useForward::Float64
end

function putParams!(algo::BlockAlgo,row::Dict)
    row[:model] = "Finkelstein PF"
    row[:mod] = "finkel"
    row[:M] = algo.M
    row[:MEquiv] = algo.M
    row[:sampType] = string(algo.sampType)
    row[:mhType] = string(algo.mhType)
    row[:hpl] = algo.histPerLoc
    row[:r] = algo.r
    row[:rsub] = algo.rsub
    row[:nIter] = algo.nIter
end

function init(algo::BlockAlgo, model::KalmanFilter)
    FrankenStep(model, algo.M, algo.r)
end

function createObservations(model, steps)
end

function saveObservations(obs, fileName)
end

function loadObservations(fileName)
end

function runAlgos(model, obs, algos, saveFile)
end

function runAlgo(model, obs, algo, saveFile, cols)
end

function stepAlgo(tau, tauhat, y, savefile, cols)
end

function finkelAlgos(sampTypes, mhTypes, histPerLocs, )
