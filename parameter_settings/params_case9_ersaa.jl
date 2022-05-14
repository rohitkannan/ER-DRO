# parameters for ER-SAA + Ridge with d_x = 100

# model and data generation files
const modelFile = "portfolio_10.jl"
const dataFile = "returns_model.jl"

# parameters for estimating the true optimal objective value using MC sampling
const numMCScenarios = 20000
const numMCReplicates = 30

# number of covariates
const numCovariates = 100

# regression method
const regressionMethod = "ridge"

# determine number of samples depending on the number of covariates
const numDataSamples = Int64[303,505,1010,2020]

# number of replicates of covariate realizations
const numCovariateReplicates = 20

# degree of nonlinearity in demand model
const degree = Float64[1,0.5,2]

# scaling factor for the additive errors
const common_errors_scaling = 0.02
const individual_errors_scaling = 0.025

# starting seed
const startingSeed = 5079
