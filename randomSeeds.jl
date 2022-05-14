# generate random seeds for use within the method

using Distributions

# random seeds for Monte Carlo estimation of optimality gap
srand(2176)
const randomSeeds_MC = ceil.(Int64,rand(Uniform(0,10000),maxNumMCReplicates))

# random seeds for constructing different model replicates
srand(4007)
const randomSeeds_models = ceil.(Int64,rand(Uniform(0,10000),maxNumModelReplicates))

# random seeds for constructing different data replicates for a given model
srand(1845)
const randomSeeds_data = ceil.(Int64,rand(Uniform(0,10000),maxNumDataReplicates))

# random seeds for constructing different covariate replicates for a given data replicate
srand(2769)
const randomSeeds_covariate = ceil.(Int64,rand(Uniform(0,10000),maxNumCovariateReplicates))
