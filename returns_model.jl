# Generate data for the random model returns

# function to generate random correlation matrices. based on https://stats.stackexchange.com/questions/124538/how-to-generate-a-large-full-rank-random-correlation-matrix-with-some-strong-cor
function generateRandomCorrMat(dim::Int64)

	betaparam::Float64 = 2.0

	partCorr = zeros(Float64,dim,dim)
	corrMat = eye(dim)

    for k = 1:dim-1
        for i = k+1:dim
			partCorr[k,i] = ((rand(Beta(betaparam,betaparam),1))[1] - 0.5)*2.0
            p::Float64 = partCorr[k,i]
            for j = (k-1):-1:1
                p = p*sqrt((1-partCorr[j,i]^2)*(1-partCorr[j,k]^2)) + partCorr[j,i]*partCorr[j,k]
            end
			corrMat[k,i] = p
			corrMat[i,k] = p
        end
    end

    permut = randperm(dim)
    corrMat = corrMat[permut, permut]

    return corrMat
end


# generate base parameters for the returns model
function generateBaseReturnsModel(numCovariates::Int64,degNum::Int64)

	# generate a common set of coefficients across all replicates
	covariate_mean = zeros(Float64, maxNumCovariates)
	covariate_covMat = generateRandomCorrMat(maxNumCovariates) + 1E-09*eye(maxNumCovariates)

	scaling_factor::Float64 = 1.0
	if degNum == 1
		scaling_factor = 1.0/0.798
	elseif degNum == 2
		scaling_factor = 1.0/0.822
	end

	# generate coefficients for the returns model
	alpha = [0.005*i for i=1:numAssets]
	beta = ((scaling_factor*[0.0125,0.0075,0.005])*(1:1:numAssets)')' 

	pert_frac::Float64 = 0.2
	for i = 1:numAssets
		for j = 2:3
			beta[i,j] *= (1 - pert_frac + 2*pert_frac*rand())
		end
		beta[i,1] = scaling_factor*0.025*i - beta[i,2] - beta[i,3]
	end
	
	if(numCovariates < 3)
		throw(ErrorException("Expected numCovariates >= 3!"))
	end

	coeff_true = zeros(Float64,numCovariates+1,numAssets)
	coeff_true[1,:] = alpha
	coeff_true[2,:] = beta[:,1]
	coeff_true[3,:] = beta[:,2]
	coeff_true[4,:] = beta[:,3]
	
	return covariate_mean, covariate_covMat, coeff_true
end



# generate returns and covariate data based on a sparse (non)linear model
function generateReturnsData(numCovariates::Int64,numSamples::Int64,degree::Float64,covariate_mean::Array{Float64},covariate_covMat::Array{Float64,2},coeff_true::Array{Float64,2})

	covariate_mean_case = covariate_mean[1:numCovariates]
	covariate_covMat_case = covariate_covMat[1:numCovariates,1:numCovariates]
	covariate_covMat_chol = (cholfact(Hermitian(covariate_covMat_case)))[:L]

	#*===========================================
	# first, construct the the covariates	
	# first covariate is simply the constant one (for the intercept)
	# the next set of covariates are iid samples from the data generation process
	
	covariate_data = zeros(Float64, numSamples, numCovariates+1)

	covariate_data[:,1] = ones(Float64,numSamples)
	covariate_tmp = rand(Normal(),maxNumCovariates,maxNumDataSamples)
	for s = 1:numSamples
		covariate_data[s,2:numCovariates+1] = covariate_covMat_chol*covariate_tmp[1:numCovariates,s] + covariate_mean_case
	end
	covariate_data = abs.(covariate_data)

	
	#*===========================================
	# next, construct the random returns from the covariates using a (sparse) (non)linear model
	
	return_errors = individual_errors_scaling*rand(Normal(), maxNumDataSamples, numAssets).*((repeat(1:1:numAssets,outer=[1,maxNumDataSamples]))')
	common_return_errors = common_errors_scaling*rand(Normal(), maxNumDataSamples)*ones(1,numAssets)
	return_data = ((covariate_data).^(degree))*coeff_true + return_errors[1:numSamples,:] + common_return_errors[1:numSamples,:]
	
	
	return return_data, covariate_data
end



# generate new realizations of the covariates
function generateCovariateReal(numCovariates::Int64,covariate_mean::Array{Float64},covariate_covMat::Array{Float64,2})

	covariate_mean_case = covariate_mean[1:numCovariates]
	covariate_covMat_case = covariate_covMat[1:numCovariates,1:numCovariates]
	covariate_covMat_chol = (cholfact(Hermitian(covariate_covMat_case)))[:L]
	
	covariate_tmp = rand(Normal(),maxNumCovariates)
	covariate_obs = ones(Float64, numCovariates+1)
	covariate_obs[2:numCovariates+1] = covariate_covMat_chol*(covariate_tmp[1:numCovariates]) + covariate_mean_case
	covariate_obs = abs.(covariate_obs)
	
	return covariate_obs
end



# generate samples from the true conditional distribution
# assume that the first entry of the covariate vector is simply the number one
function generateTrueCondScenarios(numScenarios::Int64,covariate_obs::Array{Float64,1},degree::Float64,coeff_true::Array{Float64,2})

	errors1 = rand(Normal(), maxNumMCScenarios, numAssets)
	errors2 = rand(Normal(), maxNumMCScenarios)

	return_errors = individual_errors_scaling*errors1[1:numScenarios,:].*((repeat(1:1:numAssets,outer=[1,numScenarios]))')
	common_return_errors = common_errors_scaling*errors2[1:numScenarios]*ones(1,numAssets)
	return_scen = repeat((((covariate_obs).^(degree))'*coeff_true)', outer = [1,numScenarios])' + return_errors + common_return_errors
	
	return return_scen
end
