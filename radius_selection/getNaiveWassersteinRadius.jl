
# Determine a covariate-independent Wasserstein radius (Algorithm 1) using a naive SAA-based DRO problem
function getWassersteinRadius(returns::Array{Float64,2},wassersteinRadii::Array{Float64,1})

	numSamples::Int64 = size(returns,1)	
	numRadii::Int64 = size(wassersteinRadii,1)
	numFolds::Int64 = 5
	
	# randomly permute the returns data
	rand_perm = randperm(numSamples)
	returns_perm = zeros(Float64,numSamples,numAssets)
	for ass = 1:numAssets
		returns_perm[:,ass] = returns[rand_perm,ass]
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
	

	# solve the naive SAA-based DRO problem with the data for each fold omitted
	solnCostsRadii = zeros(Float64,numRadii)
	for f = 1:numFolds
		scenarios_fold = union([1:index_folds[f,1]-1,index_folds[f,2]+1:numSamples]...)
		returns_scen_fold = returns_perm[scenarios_fold,:]
		for r = 1:numRadii
			z_soln_wass, ~ = solveWassersteinModel(returns_scen_fold,wassersteinRadii[r]) 
			solnCostsRadii[r] += estimateCostOfSoln(z_soln_wass,returns_perm[index_folds[f,1]:index_folds[f,2],:])/numFolds
		end
	end

	return wassersteinRadii[indmin(solnCostsRadii)], solnCostsRadii
end

