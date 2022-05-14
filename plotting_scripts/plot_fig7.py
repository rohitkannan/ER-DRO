# Python script for plotting the radii of the different Wasserstein ER-DRO + OLS algorithms for d_x = 100
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from statistics import mean, quantiles
from scipy import stats


num_degrees = 3
degrees = [1.0, 0.5, 2.0]
num_covariate_dimensions = 1
covariate_dimensions = np.array([100], dtype=np.int32)
num_data_replicates = 50
num_covariate_replicates = 20
num_saa_replicates = 30
num_sample_sizes = 4
sample_sizes = np.array([[505, 1010, 2020, 5050]], dtype=np.int32)
num_risk_levels = 28
conf_level = 0.99
t_value = stats.t.ppf(conf_level, num_saa_replicates-1)
wass_radii = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
whiskers = (2.0,98.0)
font_size = 32
line_width = 4
fig_size = (10,8)
dpi_val = 1200

# First load the FI-SAA data with 100K scenarios
SAA_data = np.zeros((num_degrees, num_covariate_replicates, num_saa_replicates))

for deg in range(num_degrees):
	for cov_rep in range(num_covariate_replicates):
		SAA_file = open("case1_saa/mod_1/deg_" + str(deg+1) + "/cov_" + str(cov_rep+1) + "/saa_obj.txt", "r")
		SAA_list = SAA_file.read().splitlines()
		SAA_data[deg,cov_rep,:] = [float(i) for i in SAA_list]
		SAA_file.close()
		
SAA_mean = np.mean(SAA_data, axis=2)



# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 3 is used to specify the radius
wass_ERDRO_doubleest_radii = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[1, 2, 3, 4]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			missing_data = []
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_doubleest/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_radius.txt") as wass_file:
					wass_list = wass_file.read().splitlines() 
					if len(wass_list) != num_covariate_replicates:
						missing_data.append(data_rep)
						print(data_rep)
					else:
						wass_ERDRO_doubleest_radii[cov_dim,deg,:,data_rep,samp_size] = [float(i) for i in wass_list]
				wass_file.close()
			if missing_data:
				exist_data = set(range(num_data_replicates))
				exist_data = list(exist_data.difference(set(missing_data)))
				for data_rep in missing_data:
					for cov_rep in range(num_covariate_replicates):
						wass_ERDRO_doubleest_radii[cov_dim,deg,cov_rep,data_rep,samp_size] = np.median([wass_ERDRO_doubleest_radii[cov_dim,deg,cov_rep,ind,samp_size] for ind in exist_data])


# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 2 is used to specify the radius
wass_ERDRO_scrambled_radii = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[1, 2, 3, 4]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_scrambled/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_radius.txt") as wass_file:
					wass_list = wass_file.read()
					wass_ERDRO_scrambled_radii[cov_dim,deg,data_rep,samp_size] = float(wass_list)
				wass_file.close()


# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 1 is used to specify the radius
wass_ERDRO_naive_radii = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[0, 1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_naive/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_radius.txt") as wass_file:
					wass_list = wass_file.read()
					wass_ERDRO_naive_radii[cov_dim,deg,data_rep,samp_size] = float(wass_list)
				wass_file.close()


# Next, load the Wasserstein ER-DRO + OLS data for different choices of the radius
wass_ERDRO_tailored_data = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_risk_levels, num_saa_replicates))
wass_ERDRO_tailored_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_risk_levels))
wass_ERDRO_tailored_ucb_dep = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_ERDRO_tailored_indep_crit = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes, num_risk_levels))
wass_ERDRO_tailored_ucb_indep = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_ERDRO_tailored_indep_radii = np.zeros((num_covariate_dimensions, num_degrees, num_data_replicates, num_sample_sizes))
wass_ERDRO_tailored_dep_radii = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[1, 2, 3, 4]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				for risk_num in range(num_risk_levels):
					with open("case" + str(3+cov_dim) + "_wass_tailored/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/risk_" + str(risk_num+1) + "/wass_obj.txt") as wass_file:
						wass_list = wass_file.read().splitlines() 
						for cov_rep in range(num_covariate_replicates):
							wass_ERDRO_tailored_data[cov_dim,deg,cov_rep,data_rep,samp_size,risk_num,:] = [float(i) for i in wass_list[cov_rep*num_saa_replicates:(cov_rep+1)*num_saa_replicates]]
							wass_ERDRO_tailored_gap = wass_ERDRO_tailored_data[cov_dim,deg,cov_rep,data_rep,samp_size,risk_num,:] - SAA_data[deg,cov_rep,:]
							wass_ERDRO_tailored_gap_mean = np.mean(wass_ERDRO_tailored_gap)
							wass_ERDRO_tailored_gap_var = np.var(wass_ERDRO_tailored_gap)
							wass_ERDRO_tailored_ucb[cov_dim,deg,cov_rep,data_rep,samp_size,risk_num] = (100.0)*(wass_ERDRO_tailored_gap_mean + t_value*np.sqrt(wass_ERDRO_tailored_gap_var/num_saa_replicates))
					wass_file.close()
					wass_ERDRO_tailored_indep_crit[cov_dim,deg,data_rep,samp_size,risk_num] = np.median(wass_ERDRO_tailored_ucb[cov_dim,deg,:,data_rep,samp_size,risk_num])
				opt_risk_num = np.argmin(wass_ERDRO_tailored_indep_crit[cov_dim,deg,data_rep,samp_size,:])
				wass_ERDRO_tailored_indep_radii[cov_dim,deg,data_rep,samp_size] = wass_radii[opt_risk_num]
				wass_ERDRO_tailored_ucb_indep[cov_dim,deg,:,data_rep,samp_size] = wass_ERDRO_tailored_ucb[cov_dim,deg,:,data_rep,samp_size,opt_risk_num]
				for cov_rep in range(num_covariate_replicates):
					opt_risk_num = np.argmin(wass_ERDRO_tailored_ucb[cov_dim,deg,cov_rep,data_rep,samp_size,:])
					wass_ERDRO_tailored_dep_radii[cov_dim,deg,cov_rep,data_rep,samp_size] = wass_radii[opt_risk_num]
					wass_ERDRO_tailored_ucb_dep[cov_dim,deg,cov_rep,data_rep,samp_size] = wass_ERDRO_tailored_ucb[cov_dim,deg,cov_rep,data_rep,samp_size,opt_risk_num]


# Plot the radii of the different Wasserstein ER-DRO + OLS algorithms
box_colors = ['magenta','black','limegreen','red','blue']
plot_y_max = [0.25, 0.2, 0.61]
for deg in range(num_degrees):
	for cov in range(num_covariate_dimensions):
		plt.figure(figsize=fig_size, dpi=dpi_val)
		plt.rcParams["font.weight"] = "bold"
		plt.rcParams["axes.labelweight"] = "bold"
		plt.rcParams["axes.titleweight"] = "bold"
		plt.rcParams.update({'font.size': font_size})

		for meth in range(5):
			data = np.zeros((num_data_replicates*num_covariate_replicates,num_sample_sizes))
			box_positions = np.zeros(num_sample_sizes)
			box_width = 0.5
			for samp_size in range(num_sample_sizes):
				if meth == 0:
					data[:,samp_size] = wass_ERDRO_tailored_dep_radii[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1 + 2.7*samp_size
				elif meth == 1:
					data[:,samp_size] = wass_ERDRO_doubleest_radii[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1.5 + 2.7*samp_size
				elif meth == 2:
					data[:,samp_size] = np.repeat(wass_ERDRO_tailored_indep_radii[cov,deg,:,samp_size].flatten(),num_covariate_replicates)
					box_positions[samp_size] = 2 + 2.7*samp_size
				elif meth == 3:
					data[:,samp_size] = np.repeat(wass_ERDRO_naive_radii[cov,deg,:,samp_size].flatten(),num_covariate_replicates)
					box_positions[samp_size] = 2.5 + 2.7*samp_size
				elif meth == 4:
					data[:,samp_size] = np.repeat(wass_ERDRO_scrambled_radii[cov,deg,:,samp_size].flatten(),num_covariate_replicates)
					box_positions[samp_size] = 3 + 2.7*samp_size

			plt.boxplot(data,whis=whiskers, showfliers=False, 
						positions=box_positions, widths=box_width,
						boxprops={'linewidth':line_width,'color':box_colors[meth]}, 
						medianprops={'linewidth':line_width,'color':box_colors[meth]},
						whiskerprops={'linewidth':line_width,'color':box_colors[meth]},
						capprops={'linewidth':line_width,'color':box_colors[meth]})
			if deg == 0:
				plt.ylabel('Wasserstein radius')
			plt.title(r'$\theta$' + ' = ' + str(degrees[deg]))
			x1,x2,y1,y2 = plt.axis()  
			plt.axis((x1,2.7*(num_sample_sizes)+0.7,0,plot_y_max[deg]))


			box_positions_overall = np.zeros(6*num_sample_sizes)
			box_labels_overall = np.empty(6*num_sample_sizes,dtype=object)
			for samp_size in range(num_sample_sizes):
				box_positions_overall[6*samp_size] = 1 + 2.7*samp_size
				box_labels_overall[6*samp_size] = r'D$\!^*$'
				box_positions_overall[6*samp_size+1] = 1.5 + 2.7*samp_size
				box_labels_overall[6*samp_size+1] = '3'
				box_positions_overall[6*samp_size+2] = 1.99 + 2.7*samp_size
				box_labels_overall[6*samp_size+2] = '\nn=' + str(sample_sizes[cov][samp_size])
				box_positions_overall[6*samp_size+3] = 2 + 2.7*samp_size
				box_labels_overall[6*samp_size+3] = r'I$\!^*$'
				box_positions_overall[6*samp_size+4] = 2.5 + 2.7*samp_size
				box_labels_overall[6*samp_size+4] = '1'
				box_positions_overall[6*samp_size+5] = 3 + 2.7*samp_size
				box_labels_overall[6*samp_size+5] = '2'
			
			plt.xticks(box_positions_overall, box_labels_overall)
			plt.tick_params(axis='x',
							which='both',
							bottom=False)

			if deg == 0:
				plt.subplots_adjust(left=0.17,right=0.99,bottom=0.12)
			else:
				plt.subplots_adjust(left=0.13,right=0.99,bottom=0.12)
			plt.savefig("figures/fig7_cov" + str(cov+3) + "_deg" + str(deg+1) + ".eps", format='eps')
