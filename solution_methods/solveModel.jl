# given scenarios of the returns, solve the SAA problem
function solveSAAModel(returns::Array{Float64,2})

	numScenarios::Int64 = size(returns,1)
	mean_returns = mean(returns,1)
	
	mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,Threads=maxNumThreads))

	@variable(mod, z[1:numAssets] >= 0)
	@variable(mod, theta[1:numScenarios] >= 0)
	@variable(mod, gamma)

	
	@objective(mod, Min, cvar_mult_factor*gamma - sum(mean_returns[i]*z[i] for i=1:numAssets) + (cvar_mult_factor/((1.0-cvar_risk_param)*numScenarios))*sum(theta[scen] for scen=1:numScenarios))
	
	@constraint(mod, sum(z[i] for i=1:numAssets) == 1)
	@constraint(mod, [scen=1:numScenarios], theta[scen] >= -sum(returns[scen,i]*z[i] for i=1:numAssets) - gamma)

	
	status = solve(mod)
	z_soln = getvalue(z)
	objValue_soln::Float64 = getobjectivevalue(mod)
	
	return z_soln, objValue_soln
end
