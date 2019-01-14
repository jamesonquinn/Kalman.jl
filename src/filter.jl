#######
# Universal Kalman filtering methods
import StatsBase.predict
import StatsBase.predict!

function predict!(kf::KalmanFilter)
    kf.x = ap(kf.f,kf.x)
    kf
end

function predict(kf::KalmanFilter)
    predict!(copy(kf))
end

function update(kf::KalmanFilter,y::Observation)
    doUpdate!(copy(kf),y)
end

function doUpdate!(kf::KalmanFilter,y::Observation)
    #println()
    (res,ph,s) = covs(kf,y)





    xn = kf.x.x + ph * (s\res)
    pn = kf.x.p - ph * (s'\ph')

    # This is an ugly hack which works for now
    if typeof(kf.x) <: AbstractUnscentedState
        kf.x = UnscentedState(xn,pn,kf.x.α,kf.x.β,kf.x.κ)
    else
        kf.x = State(xn,pn)
    end
    kf
end

function predictUpdate(kf::KalmanFilter,y::Observation)
    update(predict(kf),y)
end

function predictdoUpdate!(kf::KalmanFilter,y::Observation)
    doUpdate!(predict!(kf),y)
end
