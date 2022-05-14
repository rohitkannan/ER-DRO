# main file for solving the Hellinger ER-DRO model using Algorithm 2 for specifying the radius

using Distributions, StatsFuns, StatsBase, JuMP, Gurobi, GLMNet


const gurobi_env = Gurobi.Env()



modRepNum = parse(ARGS[1])
degNum = parse(ARGS[2])
dataRepNum = parse(ARGS[3])
sampSizeNum = parse(ARGS[4])



#*========= MODEL INFORMATION ==========
const caseNum = 1
const paramsFile = "params_case" * string(caseNum) * "_hellinger_scrambled.jl"
#*======================================



#*========= INCLUDE FILES ====================
include("maxParameters.jl")
include("randomSeeds.jl")

include(paramsFile)
include(modelFile)
include(dataFile)

include("regressionMethods.jl")
include("genERSAAScenarios.jl")
include("solveHellingerModel.jl")
include("estimateSolnCost.jl")
include("evaluateSoln.jl")
include("getScrambledHellingerRadius.jl")
#*======================================



# set seed for reproducibility
srand(startingSeed)



# directory name for storing results
const baseDirName = "case" * string(caseNum) * "_hellinger_scrambled/" * "mod_" * string(modRepNum) * "/" * "deg_" * string(degNum) * "/"
const subDirName = baseDirName * "rep_" * string(dataRepNum) * "/" * regressionMethod * "/" * "samp_" * string(sampSizeNum) * "/"
mkpath(subDirName)



#*========= PRINT OPTIONS ====================
const storeResults = true

const infoFile = "ddsp.txt"
const modelDataFile = "model_data.txt"
const numSampFile = "num_samples_" * regressionMethod * ".txt"

const chisqObjFile = "hellinger_obj.txt"
const chisqDDObjFile = "hellinger_ddobj.txt"
const chisqTimeFile = "hellinger_time.txt"
const radiusTimeFile = "radius_time.txt"
const chisqRadiusFile = "hellinger_radius.txt"
const radiusQualityFile = "radius_quality.txt"

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
		write(f,"Hellinger radius: $hellingerRadius \n")
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



#*========= GET HELLINGER RADIUS ==============

tic()

optChisqRadius, radiusQuality = getHellingerRadius(returns_data,covariate_data,hellingerRadius,regressionMethod,coeff_true)

radiusTime = toq()



#*========= STORE RESULTS ====================
if(storeResults) 

	time_radius_file = subDirName * radiusTimeFile
	open(time_radius_file, "w") do f
		write(f,"$radiusTime") 
	end

	radius_file = subDirName * chisqRadiusFile		
	open(radius_file, "w") do f 
		write(f,"$optChisqRadius")
	end

	radius_quality_file = subDirName * radiusQualityFile
	open(radius_quality_file, "w") do f
		for i = 1:length(radiusQuality)
			write(f,"$(radiusQuality[i]) \n")
		end
	end
end
#*============================================



#*========= GENERATE COVARIATE REPLICATE ============== 


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



#*========= SOLVE HELLINGER ER-DRO MODEL ====================

	tic()


	returns_scen_ersaa, ~ = generateERSAAScenarios(returns_data,covariate_data,covariate_obs,regressionMethod,coeff_true)

	z_soln_chisq, objDDChisq = solveHellingerModel(returns_scen_ersaa,optChisqRadius)

	chisqObjEstimates = estimateSolnQuality(z_soln_chisq,covariate_obs,degree[degNum],numMCScenarios,numMCReplicates,coeff_true)


	chisqTime = toq()


#*========= STORE RESULTS ====================
	if(storeResults)
	
		# write data to text file
		obj_est_file = subDirName * chisqObjFile
		open(obj_est_file, "a") do f
			for i = 1:length(chisqObjEstimates)
				write(f,"$(chisqObjEstimates[i]) \n")
			end
		end	
	
		ddobj_est_file = subDirName * chisqDDObjFile
		open(ddobj_est_file, "a") do f
			write(f,"$objDDChisq \n")
		end	
	
		time_obj_file = subDirName * chisqTimeFile
		open(time_obj_file, "a") do f
			write(f,"$chisqTime \n")
		end	

	end
#*============================================


end # covRepNum
