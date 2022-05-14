# main file for solving the J-SAA and J+-SAA problems

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi, GLMNet


const gurobi_env = Gurobi.Env()



modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
sampSizeNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 3
const paramsFile = "params_case" * string(caseNum) * "_jsaa.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
include(modelFile)
include(dataFile)

include("regressionMethods.jl")
include("genJSAAScenarios.jl")
include("solveModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_jsaa/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/" * regressionMethod * "/" * "samp_" * string(sampSizeNum) * "/"
mkpath(subDirName)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples_" * regressionMethod * ".txt"

const jsaaObjFile = "jsaa_obj.txt"
const jsaaDDObjFile = "jsaa_ddobj.txt"
const jpsaaObjFile = "jpsaa_obj.txt"
const jpsaaDDObjFile = "jpsaa_ddobj.txt"
const jsaaTimeFile = "jsaa_time.txt"

const covRealFile = "covariate_obs.txt"
#*======================================
	


#*========= GENERATE MODEL PARAMETERS ====================
	
srand(randomSeeds_models[modRepNum])

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
		write(f,"sample sizes: $numDataSamples \n")
		write(f,"common_errors_scaling: $common_errors_scaling \n")
		write(f,"individual_errors_scaling: $individual_errors_scaling \n")
		write(f,"regressionMethod: $regressionMethod \n\n")
		
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
	
	samp_size_file = baseDirName * numSampFile
	open(samp_size_file, "w") do f
		for i = 1:length(numDataSamples)
			write(f,"$(numDataSamples[i]) \n")
		end
	end
	
end
#*============================================



#*========= GENERATE DATA REPLICATE ====================

srand(randomSeeds_data[dataRepNum])

returns_data, covariate_data = generateReturnsData(numCovariates,numDataSamples[sampSizeNum],degree[degNum],covariate_mean,covariate_covMat,coeff_true)



#*========= GENERATE COVARIATE REPLICATE ====================

for covRepNum = 1:numCovariateReplicates


	srand(randomSeeds_covariate[covRepNum])

	covariate_obs = generateCovariateReal(numCovariates,covariate_mean,covariate_covMat)



#*========= STORE RESULTS ====================
	if(storeResults)
	
		cov_real_file = baseDirName * covRealFile
		open(cov_real_file, "a") do f
			write(f,"covariate_obs = $covariate_obs \n")
		end	
	
	end
#*============================================



#*========= SOLVE JSAA MODEL ====================

	tic()


	returns_scen_jsaa, returns_scen_jpsaa = generateJSAAScenarios(returns_data,covariate_data,covariate_obs,regressionMethod)

	z_soln_jsaa, objDDJSAA = solveSAAModel(returns_scen_jsaa)

	jsaaObjEstimates = estimateSolnQuality(z_soln_jsaa,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


        z_soln_jpsaa, objDDJpSAA = solveSAAModel(returns_scen_jpsaa)

        jpsaaObjEstimates = estimateSolnQuality(z_soln_jpsaa,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


	jsaaTime = toq()


#*========= STORE RESULTS ====================
	if(storeResults)
	
		# write data to text file
		obj_est_file1 = subDirName * jsaaObjFile
		open(obj_est_file1, "a") do f
			for i = 1:length(jsaaObjEstimates)
				write(f,"$(jsaaObjEstimates[i]) \n")
			end
		end

		obj_est_file2 = subDirName * jpsaaObjFile
		open(obj_est_file2, "a") do f
			for i = 1:length(jpsaaObjEstimates)
				write(f,"$(jpsaaObjEstimates[i]) \n")
			end
		end
	
		ddobj_est_file1 = subDirName * jsaaDDObjFile
		open(ddobj_est_file1, "a") do f
			write(f,"$objDDJSAA \n")
		end

		ddobj_est_file2 = subDirName * jpsaaDDObjFile
		open(ddobj_est_file2, "a") do f
			write(f,"$objDDJpSAA \n")
		end
	
		time_obj_file = subDirName * jsaaTimeFile
		open(time_obj_file, "a") do f
			write(f,"$jsaaTime \n")
		end	
	
	end
#*============================================


end # covariate replicates
