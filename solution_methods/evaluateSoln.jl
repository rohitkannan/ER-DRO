# ROUTINES FOR EVALUATING QUALITY OF COMPUTED FIRST-STAGE SOLUTION


# estimate the objective value of a given first-stage solution
function estimateSolnQuality(z_soln::Array{Float64},covariate_obs::Array{Float64},degree::Float64,numMCScenarios::Int64,numMCReplicates::Int64,coeff_true::Array{Float64,2})
	
	objEstimates = zeros(Float64,numMCReplicates)

	for mc = 1:numMCReplicates
		srand(randomSeeds_MC[mc])
	
		returns_scen_MC = generateTrueCondScenarios(numMCScenarios,covariate_obs,degree,coeff_true)
		
		# estimate the out-of-sample cost of given solution
		objEstimates[mc] = estimateCostOfSoln(z_soln,returns_scen_MC)
	end
	
	return objEstimates
end
