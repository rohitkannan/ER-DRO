# parameters for the Full Information SAA (FI-SAA) problem with 100K scenarios

# model and data generation files
const modelFile = "portfolio_10.jl"
const dataFile = "returns_model.jl"

# parameters for estimating the true optimal objective value using MC sampling
const numMCScenarios = 100000
const numMCReplicates = 30

# degree of nonlinearity in demand model
const degree = Float64[1,0.5,2]

# scaling factor for the additive errors
const common_errors_scaling = 0.02
const individual_errors_scaling = 0.025

# starting seed
const startingSeed = 5079
