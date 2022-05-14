# Python script for plotting ER-SAA + OLS and Wasserstein ER-DRO + OLS (using Algorithm 2) in-sample objective values
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
whiskers = (5.0,95.0)
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


# Next, load the ER-SAA + OLS data
ERSAA_cert = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([[1, 2, 3, 4]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(ERSAA_samp_num[cov_dim][samp_size]+1) + "/ersaa_ddobj.txt") as ERSAA_file:
					ERSAA_list = ERSAA_file.read().splitlines() 
					ERSAA_cert[cov_dim,deg,:,data_rep,samp_size] = [(100.0)*(float(ERSAA_list[i]) - SAA_mean[deg,i]) for i in range(len(ERSAA_list))]
				ERSAA_file.close()


# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 2 is used to specify the radius
wass_ERDRO_scrambled_cert = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[1, 2, 3, 4]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_scrambled/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_ddobj.txt") as wass_file:
					wass_list = wass_file.read().splitlines() 
					wass_ERDRO_scrambled_cert[cov_dim,deg,:,data_rep,samp_size] = [(100.0)*(float(wass_list[i]) - SAA_mean[deg,i]) for i in range(len(wass_list))]
				wass_file.close()

# Plot ER-SAA + OLS and Wasserstein ER-DRO + OLS (using Algorithm 2) in-sample objective values
box_colors = ['blue', 'red']
plot_y_min = [-200, -200, -250]
plot_y_max = [200, 125, 500]
for deg in range(num_degrees):
	for cov in range(num_covariate_dimensions):
		plt.figure(figsize=fig_size, dpi=dpi_val)
		plt.rcParams["font.weight"] = "bold"
		plt.rcParams["axes.labelweight"] = "bold"
		plt.rcParams["axes.titleweight"] = "bold"
		plt.rcParams.update({'font.size': font_size})

		for meth in range(2):
			data = np.zeros((num_data_replicates*num_covariate_replicates,num_sample_sizes))
			box_positions = np.zeros(num_sample_sizes)
			box_width = 0.5
			for samp_size in range(num_sample_sizes):
				if meth == 0:
					data[:,samp_size] = wass_ERDRO_scrambled_cert[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1 + 1.2*samp_size
				elif meth == 1:
					data[:,samp_size] = ERSAA_cert[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1.5 + 1.2*samp_size

			plt.boxplot(data,whis=whiskers, showfliers=False, 
						positions=box_positions, widths=box_width,
						boxprops={'linewidth':line_width,'color':box_colors[meth]}, 
						medianprops={'linewidth':line_width,'color':box_colors[meth]},
						whiskerprops={'linewidth':line_width,'color':box_colors[meth]},
						capprops={'linewidth':line_width,'color':box_colors[meth]})
			if deg == 0:
				plt.ylabel('Normalized certificate %')
			plt.title(r'$\theta$' + ' = ' + str(degrees[deg]))
			x1,x2,y1,y2 = plt.axis()  
			plt.axis((x1,1.2*(num_sample_sizes)+0.7,plot_y_min[deg],plot_y_max[deg]))


			box_positions_overall = np.zeros(3*num_sample_sizes)
			box_labels_overall = np.empty(3*num_sample_sizes,dtype=object)
			for samp_size in range(num_sample_sizes):
				box_positions_overall[3*samp_size] = 1 + 1.2*samp_size
				box_labels_overall[3*samp_size] = 'W'
				if samp_size == 0:
					box_positions_overall[3*samp_size+1] = 1.05 + 1.2*samp_size
				elif samp_size == 1:
					box_positions_overall[3*samp_size+1] = 1.1 + 1.2*samp_size
				elif samp_size == 2:
					box_positions_overall[3*samp_size+1] = 1.25 + 1.2*samp_size
				elif samp_size == 3:
					box_positions_overall[3*samp_size+1] = 1.35 + 1.2*samp_size
				box_labels_overall[3*samp_size+1] = '\nn=' + str(sample_sizes[cov][samp_size])
				box_positions_overall[3*samp_size+2] = 1.5 + 1.2*samp_size
				box_labels_overall[3*samp_size+2] = 'E'
			
			plt.xticks(box_positions_overall, box_labels_overall)
			plt.tick_params(axis='x',
							which='both',
							bottom=False)

			if deg == 0:
				plt.subplots_adjust(left=0.2,right=0.99,bottom=0.12)
			else:
				plt.subplots_adjust(left=0.14,right=0.99,bottom=0.12)
			plt.savefig("figures/fig4_cov" + str(cov+3) + "_deg" + str(deg+1) + ".eps", format='eps')
