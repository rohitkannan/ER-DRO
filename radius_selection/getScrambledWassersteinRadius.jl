
# Determine a covariate-independent Wasserstein radius (Algorithm 2) using a ER-SAA-based DRO problem
function getWassersteinRadius(returns_data::Array{Float64,2},covariates_data::Array{Float64,2},wassersteinRadii::Array{Float64,1},regressionMethod::String,coeff_true::Array{Float64,2})

	numSamples::Int64 = size(returns_data,1)	
	numRadii::Int64 = size(wassersteinRadii,1)
	numFolds::Int64 = 5
	numCovariatePointsForRadius::Int64 = 50
	
	# randomly permute the returns and covariate data
	rand_perm = randperm(numSamples)
	returns_perm = zeros(Float64,numSamples,numAssets)
	covariates_perm = zeros(Float64,numSamples,numCovariates+1)
	for ass = 1:numAssets
		returns_perm[:,ass] = returns_data[rand_perm,ass]
	end
	for cov = 1:numCovariates+1
		covariates_perm[:,cov] = covariates_data[rand_perm,cov]
	end

	
	# split the data into numFolds folds
	numBaseElements::Int64 = floor.(Int64,numSamples*1.0/numFolds)
	numRemElements::Int64 = numSamples - numFolds*numBaseElements
	index_folds = zeros(Int64,numFolds,2)
	index_folds[1,1] = 1
	index_folds[1,2] = numBaseElements
	if(numRemElements > 0)
		index_folds[1,2] += 1
		numRemElements -= 1
	end
	for f = 2:numFolds
		index_folds[f,1] = index_folds[f-1,2] + 1
		index_folds[f,2] = index_folds[f,1] + (numBaseElements-1)
		if(numRemElements > 0)
			index_folds[f,2] += 1
			numRemElements -= 1
		end
	end
	

	# solve the ER-SAA-based DRO problem with the data for each fold omitted
	solnCostsRadii = zeros(Float64,numRadii)
	numCovariatePoints::Int64 = min(numBaseElements,numCovariatePointsForRadius)
	for f = 1:numFolds
		scenarios_fold = union([1:index_folds[f,1]-1,index_folds[f,2]+1:numSamples]...)
		returns_data_fold = returns_perm[scenarios_fold,:]
		covariates_data_fold = covariates_perm[scenarios_fold,:]
		covariate_values_fold = covariates_perm[index_folds[f,1]:index_folds[f,2],:]
		for cov = 1:numCovariatePoints
			returns_scen_fold, ~ = generateERSAAScenarios(returns_data_fold,covariates_data_fold,covariate_values_fold[cov,:],regressionMethod,coeff_true)
		
			for r = 1:numRadii
				z_soln_wass, ~ = solveWassersteinModel(returns_scen_fold,wassersteinRadii[r]) 
				solnCostsRadii[r] += (estimateCostOfSoln(z_soln_wass,returns_perm[index_folds[f,1]:index_folds[f,2],:]))/(numFolds*numCovariatePoints)
			end			
		end
	end

	return wassersteinRadii[indmin(solnCostsRadii)], solnCostsRadii
end

