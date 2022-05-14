# Estimate cost of decision given scenarios of the returns

# given decisions and scenarios of the returns, estimate the objective cost
function estimateCostOfSoln(z_soln::Array{Float64},returns::Array{Float64,2})

	numScenarios::Int64 = size(returns,1)
	mean_returns = mean(returns,1)
	
	solnCost::Float64 = -sum(mean_returns[i]*z_soln[i] for i=1:numAssets)
	

	mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0,Threads=maxNumThreads))
	
	@variable(mod, gamma)
	@variable(mod, theta[1:numScenarios] >= 0)

	
	@objective(mod, Min, gamma + (1.0/((1.0-cvar_risk_param)*numScenarios))*sum(theta[scen] for scen=1:numScenarios))
	
	@constraint(mod, [scen=1:numScenarios], theta[scen] >= -sum(returns[scen,i]*z_soln[i] for i=1:numAssets) - gamma)
	

	status = solve(mod)
	solnCost += cvar_mult_factor*getobjectivevalue(mod)
	
	return solnCost
end
