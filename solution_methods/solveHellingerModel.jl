# Solve a SAA of the portfolio model given scenarios of the returns
include("solveModel.jl")


# given scenarios of the returns, solve the Hellinger ER-DRO problem
function solveHellingerModel(returns::Array{Float64,2},radius::Float64)

	z_soln = zeros(Float64,numAssets)
	objValue_soln::Float64 = Inf

	if(radius > 0.0)

		numScenarios::Int64 = size(returns,1)
	
		mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,Threads=maxNumThreads))

		@variable(mod, z[1:numAssets] >= 0)
		@variable(mod, h[1:numScenarios])
		@variable(mod, theta[1:numScenarios])
		@variable(mod, lambda >= 0)
		@variable(mod, mu)
		@variable(mod, tau)

	
		@objective(mod, Min, mu + (radius-1)*lambda + (1/numScenarios)*sum(theta[scen] for scen=1:numScenarios))
	
		@constraint(mod, sum(z[i] for i=1:numAssets) == 1)
		@constraint(mod, [scen=1:numScenarios], h[scen] >= cvar_mult_factor*tau - sum(returns[scen,i]*z[i] for i=1:numAssets))
		@constraint(mod, [scen=1:numScenarios], h[scen] >= cvar_mult_factor*(1.0 - (1.0/(1.0-cvar_risk_param)))*tau - (1.0 + (cvar_mult_factor/(1.0-cvar_risk_param)))*sum(returns[scen,i]*z[i] for i=1:numAssets))
		@constraint(mod, [scen=1:numScenarios], lambda + mu >= h[scen])
		@constraint(mod, [scen=1:numScenarios], norm([2*lambda, (theta[scen] - lambda + h[scen] - mu)]) <= theta[scen] + lambda - h[scen] + mu)

	
		status = solve(mod)
		z_soln = getvalue(z)
		objValue_soln = getobjectivevalue(mod)

	else

		z_soln, objValue_soln = solveSAAModel(returns)

	end
	
	return z_soln, objValue_soln
end
