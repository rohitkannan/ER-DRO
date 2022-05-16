# Codes for the paper "Residuals-Based Distributionally Robust Optimization with Covariate Information"

The included <a href = "https://github.com/rohitkannan/ER-DRO/blob/main/Residuals-based%20DRO%20with%20covariate%20information%20(R1).pdf" target="_blank">PDF</a> is the latest version of the paper.

Results in the paper were generated using Julia 0.6.4, Python 3.8.5, JuMP 0.18.5, Gurobi 8.1.0, GLMNet 0.3.0, and NumPy 1.19.2. Most of the code should readily port over to later versions of Julia.


**Folder Contents:**
* The main files for the different methods lie in the "main_files" folder
* Implementations of the different solution methods are in the "solution_methods" folder
* Implementations of the different radius selection strategies lie in the "radius_selection" folder
* Parameter settings for the different cases are in the "parameter_settings" folder
* The different instances for the cases are in the "instances" folder
* Shell scripts and submit files for HTCondor are in the "shell_scripts" folder
* Python scripts for plotting the results are in the "plotting_scripts" folder


**Case Details:**
* case1, case2, and case3: OLS regression for d_x = 3, d_x = 10, and d_x = 100, respectively
* case4, case5, and case6: Lasso regression for d_x = 3, d_x = 10, and d_x = 100, respectively
* case7, case8, and case9: Ridge regression for d_x = 3, d_x = 10, and d_x = 100, respectively
* "_Naive_ Radius" corresponds to Algorithm 1 for specifying the radius
* "_Scrambled_ Radius" corresponds to Algorithm 2 for specifying the radius
* "_Doubleest_ Radius" corresponds to Algorithm 3 for specifying the radius
* "_Tailored_ Radius" corresponds to choosing the radius with the "best out-of-sample performance"
