# given scenarios of the returns, solve the Sample Robust ER-DRO problem
function solveSampleRobustModel(returns::Array{Float64,2},radius::Float64)

	numScenarios::Int64 = size(returns,1)

	
	mod = Model(solver=GurobiSolver(gurobi_env,Presolve=0,OutputFlag=0))

	@variable(mod, z[1:numAssets] >= 0)
	@variable(mod, theta[1:numScenarios])
	@variable(mod, lambda)
	@variable(mod, tau)

	
	@objective(mod, Min, (1.0/numScenarios)*sum(theta[scen] for scen=1:numScenarios))
	
	@constraint(mod, sum(z[i] for i=1:numAssets) == 1)
	@constraint(mod, [scen=1:numScenarios], theta[scen] >= cvar_mult_factor*tau - sum(returns[scen,i]*z[i] for i=1:numAssets) + radius*lambda)
	@constraint(mod, [scen=1:numScenarios], theta[scen] >= cvar_mult_factor*(1.0 - (1.0/(1.0-cvar_risk_param)))*tau - (1.0 + (cvar_mult_factor/(1.0-cvar_risk_param)))*sum(returns[scen,i]*z[i] for i=1:numAssets) + (1.0 + (cvar_mult_factor/(1.0-cvar_risk_param)))*radius*lambda)
	@constraint(mod, [i=1:numAssets], z[i] <= lambda)

	
	status = solve(mod)
	z_soln = getvalue(z)
	objValue_soln::Float64 = getobjectivevalue(mod)
	
	return z_soln, objValue_soln
end
