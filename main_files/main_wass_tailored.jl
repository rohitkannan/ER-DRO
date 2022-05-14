# main file for solving the Wasserstein ER-DRO problem using different values of the radius

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi, GLMNet


const gurobi_env = Gurobi.Env()



modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
riskNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 1
const paramsFile = "params_case" * string(caseNum) * "_wass_tailored.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
include(modelFile)
include(dataFile)

include("regressionMethods.jl")
include("genERSAAScenarios.jl")
include("solveWassersteinModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_wass_tailored/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/" * regressionMethod * "/"
mkpath(subDirName)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples_" * regressionMethod * ".txt"

const wassObjFile = "wass_obj.txt"
const wassDDObjFile = "wass_ddobj.txt"
const wassTimeFile = "wass_time.txt"

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
		write(f,"Wasserstein radius: $wassersteinRadius \n")
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


for sampSizeNum = 1:length(numDataSamples)


	srand(randomSeeds_data[dataRepNum])

	returns_data, covariate_data = generateReturnsData(numCovariates,numDataSamples[sampSizeNum],degree[degNum],covariate_mean,covariate_covMat,coeff_true)


	subDirName2 = subDirName * "samp_" * string(sampSizeNum) * "/" * "risk_" * string(riskNum) * "/"
	mkpath(subDirName2)



#*========= GENERATE COVARIATE REPLICATE ====================

	for covRepNum = 1:numCovariateReplicates


		srand(randomSeeds_covariate[covRepNum])

		covariate_obs = generateCovariateReal(numCovariates,covariate_mean,covariate_covMat)

	    returns_scen_ersaa, ~ = generateERSAAScenarios(returns_data,covariate_data,covariate_obs,regressionMethod,coeff_true)



#*========= STORE RESULTS ====================
		if(storeResults)

			cov_real_file = baseDirName * covRealFile
			open(cov_real_file, "a") do f
					write(f,"covariate_obs = $covariate_obs \n")
			end

		end
#*============================================



#*========= SOLVE WASSERSTEIN ER-DRO MODEL ====================

	    tic()


	    z_soln_wass, objDDWass = solveWassersteinModel(returns_scen_ersaa,wassersteinRadius[riskNum])                                                                                                
	
		wassObjEstimates = estimateSolnQuality(z_soln_wass,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


	    wassTime = toq()



#*========= STORE RESULTS ====================
		if(storeResults)
	
			# write data to text file
			obj_est_file = subDirName2 * wassObjFile
			open(obj_est_file, "a") do f
				for i = 1:numMCReplicates
					write(f,"$(wassObjEstimates[i]) \n")
				end
			end	
	
			ddobj_est_file = subDirName2 * wassDDObjFile
			open(ddobj_est_file, "a") do f
				write(f,"$objDDWass \n")
			end	
	
			time_obj_file = subDirName2 * wassTimeFile
			open(time_obj_file, "a") do f
				write(f,"$wassTime \n")
			end	

		end
#*============================================



	end # covRepNum

end # sampSizeNum
