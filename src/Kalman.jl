module Kalman

export State,
AbstractState,
UnscentedState,
Observation, 
KalmanFilter, 
LinearKalmanFilter, 
BasicKalmanFilter, 
Model, 
ObservationModel, 
LinearModel, 
LinearObservationModel, 
predict, 
update,
predictUpdate,
predict!,
doUpdate!,
predictdoUpdate!,
BasicExtendedKalmanFilter,
NonlinearModel,
NonlinearObservationModel,
ExtendedKalmanFilter,
AdditiveUnscentedKalmanFilter,
UnscentedKalmanFilter,
AdditiveUnscentedModel,
AdditiveUnscentedObservationModel,
sigma,
AugmentedUnscentedModel,
AugmentedUnscentedObservationModel,
AugmentedUnscentedKalmanFilter,
augment

include("bkf.jl")
include("filter.jl")
include("ekf.jl")
include("unscented.jl")
include("augmented.jl")

# Nullable Observations work but only for Julia v0.4
if VERSION >= v"0.4.0-"
    include("nullable.jl")
end

end # module
