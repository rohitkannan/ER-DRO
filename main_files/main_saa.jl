# main file for solving the full-information SAA problem with a large number of scenarios

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi


const gurobi_env = Gurobi.Env()


modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
covRepNum = parse(ARGS[3])


#*========= MODEL INFORMATION ==========
const caseNum = 1
const paramsFile = "params_case" * string(caseNum) * "_saa.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
include(modelFile)
include(dataFile)

include("solveModel.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_saa/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "cov_" * string(covRepNum) * "/"
mkpath(subDirName)


#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"

const fullinfSAAObjFile = "saa_obj.txt"
const fullinfSAATimeFile = "saa_time.txt"
const covRealFile = "covariate_obs.txt"
#*======================================
	


#*========= GENERATE MODEL PARAMETERS ====================
	
srand(randomSeeds_models[modRepNum])	
numCovariates = 3

covariate_mean, covariate_covMat, coeff_true = generateBaseReturnsModel(numCovariates,degNum)


#*========= STORE RESULTS ====================
if(storeResults)
	
	# write details to text file, including some key details about the test instance
	details_file = baseDirName * infoFile
	open(details_file, "w") do f
		write(f,"case number: $caseNum \n")
		write(f,"numAssets: $numAssets \n")
		write(f,"cvar_mult_factor: $cvar_mult_factor \n")
		write(f,"cvar_risk_param: $cvar_risk_param \n")
		write(f,"model replicate number: $modRepNum \n")
		write(f,"degree: $(degree[degNum]) \n")
		write(f,"common_errors_scaling: $common_errors_scaling \n")
                write(f,"individual_errors_scaling: $individual_errors_scaling \n")
		write(f,"numMCScenarios: $numMCScenarios \n\n")
		
		write(f,"randomSeeds_MC: $randomSeeds_MC \n")
		write(f,"randomSeeds_models: $randomSeeds_models \n")
		write(f,"randomSeeds_data: $randomSeeds_data \n")
		write(f,"randomSeeds_covariate: $randomSeeds_covariate \n")
	end	
	
	mod_data_file = baseDirName * modelDataFile
	open(mod_data_file, "w") do f
		write(f,"covariate_mean = $covariate_mean \n")
		write(f,"covariate_covMat = $covariate_covMat \n")
		write(f,"coeff_true = $coeff_true \n")
	end
	
end
#*============================================



#*========= GENERATE DATA REPLICATE ====================

srand(randomSeeds_covariate[covRepNum])

covariate_obs = generateCovariateReal(numCovariates,covariate_mean,covariate_covMat)


#*========= STORE RESULTS ====================
if(storeResults)
	
	cov_real_file = subDirName * covRealFile
	open(cov_real_file, "w") do f
		write(f,"covariate_obs = $covariate_obs")
	end	

end
#*============================================



#*========= SOLVE FULL INFORMATION SAA MODEL ====================

tic()

SAAObj = zeros(Float64,numMCReplicates)

for saaRepNum = 1:numMCReplicates

	# create scenarios from the true conditional distribution
	srand(randomSeeds_MC[saaRepNum])

	returns_scen_MC = generateTrueCondScenarios(numMCScenarios,covariate_obs,degree[degNum],coeff_true)

	# solve SAA problems to estimate true objective value at this conditioning
	~, SAAObj[saaRepNum] = solveSAAModel(returns_scen_MC)

end

SAATime = toq()


#*========= STORE RESULTS ====================
if(storeResults)

	obj_est_file = subDirName * fullinfSAAObjFile
	open(obj_est_file, "w") do f
		for i = 1:numMCReplicates
			write(f,"$(SAAObj[i]) \n")
		end
	end

	time_obj_file = subDirName * fullinfSAATimeFile
	open(time_obj_file, "w") do f
		write(f,"$SAATime")
	end	
	
end
#*============================================
