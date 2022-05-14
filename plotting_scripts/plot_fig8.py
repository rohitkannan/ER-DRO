# Python script for plotting ER-SAA + OLS vs J-SAA + OLS and Wasserstein ER-DRO + OLS using Algorithms 2 and 3 for d_x = 100
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
sample_sizes = np.array([[303, 505, 1010, 2020]], dtype=np.int32)
conf_level = 0.99
t_value = stats.t.ppf(conf_level, num_saa_replicates-1)
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
ERSAA_data = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_saa_replicates))
ERSAA_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
ERSAA_samp_num = np.array([[0, 1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_ersaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(ERSAA_samp_num[cov_dim][samp_size]+1) + "/ersaa_obj.txt") as ERSAA_file:
					ERSAA_list = ERSAA_file.read().splitlines() 
					for cov_rep in range(num_covariate_replicates):
						ERSAA_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] = [float(i) for i in ERSAA_list[cov_rep*num_saa_replicates:(cov_rep+1)*num_saa_replicates]]
						ERSAA_gap = ERSAA_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] - SAA_data[deg,cov_rep,:]
						ERSAA_gap_mean = np.mean(ERSAA_gap)
						ERSAA_gap_var = np.var(ERSAA_gap)
						ERSAA_ucb[cov_dim,deg,cov_rep,data_rep,samp_size] = (100.0)*(ERSAA_gap_mean + t_value*np.sqrt(ERSAA_gap_var/num_saa_replicates))
				ERSAA_file.close()


# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 3 is used to specify the radius
wass_ERDRO_doubleest_data = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_saa_replicates))
wass_ERDRO_doubleest_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[0, 1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			missing_data = []
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_doubleest/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_obj.txt") as wass_file:
					wass_list = wass_file.read().splitlines() 
					if len(wass_list) != (num_covariate_replicates*num_saa_replicates):
						missing_data.append(data_rep)
						print(data_rep)
					else:
						for cov_rep in range(num_covariate_replicates):
							wass_ERDRO_doubleest_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] = [float(i) for i in wass_list[cov_rep*num_saa_replicates:(cov_rep+1)*num_saa_replicates]]
							wass_ERDRO_doubleest_gap = wass_ERDRO_doubleest_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] - SAA_data[deg,cov_rep,:]
							wass_ERDRO_doubleest_gap_mean = np.mean(wass_ERDRO_doubleest_gap)
							wass_ERDRO_doubleest_gap_var = np.var(wass_ERDRO_doubleest_gap)
							wass_ERDRO_doubleest_ucb[cov_dim,deg,cov_rep,data_rep,samp_size] = (100.0)*(wass_ERDRO_doubleest_gap_mean + t_value*np.sqrt(wass_ERDRO_doubleest_gap_var/num_saa_replicates))
				wass_file.close()
			if missing_data:
				exist_data = set(range(num_data_replicates))
				exist_data = list(exist_data.difference(set(missing_data)))
				for data_rep in missing_data:
					for cov_rep in range(num_covariate_replicates):
						wass_ERDRO_doubleest_ucb[cov_dim,deg,cov_rep,data_rep,samp_size] = np.median([wass_ERDRO_doubleest_ucb[cov_dim,deg,cov_rep,ind,samp_size] for ind in exist_data])


# Next, load the Wasserstein ER-DRO + OLS data when Algorithm 2 is used to specify the radius
wass_ERDRO_scrambled_data = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_saa_replicates))
wass_ERDRO_scrambled_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
wass_samp_num = np.array([[0, 1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_wass_scrambled/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(wass_samp_num[cov_dim][samp_size]+1) + "/wass_obj.txt") as wass_file:
					wass_list = wass_file.read().splitlines() 
					for cov_rep in range(num_covariate_replicates):
						wass_ERDRO_scrambled_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] = [float(i) for i in wass_list[cov_rep*num_saa_replicates:(cov_rep+1)*num_saa_replicates]]
						wass_ERDRO_scrambled_gap = wass_ERDRO_scrambled_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] - SAA_data[deg,cov_rep,:]
						wass_ERDRO_scrambled_gap_mean = np.mean(wass_ERDRO_scrambled_gap)
						wass_ERDRO_scrambled_gap_var = np.var(wass_ERDRO_scrambled_gap)
						wass_ERDRO_scrambled_ucb[cov_dim,deg,cov_rep,data_rep,samp_size] = (100.0)*(wass_ERDRO_scrambled_gap_mean + t_value*np.sqrt(wass_ERDRO_scrambled_gap_var/num_saa_replicates))
				wass_file.close()


# Next, load the J-SAA + OLS data
JSAA_data = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes, num_saa_replicates))
JSAA_ucb = np.zeros((num_covariate_dimensions, num_degrees, num_covariate_replicates, num_data_replicates, num_sample_sizes))
jack_samp_num = np.array([[0, 1, 2, 3]], dtype=np.int32)

for cov_dim in range(num_covariate_dimensions):
	for samp_size in range(num_sample_sizes):
		for deg in range(num_degrees):
			for data_rep in range(num_data_replicates):
				with open("case" + str(3+cov_dim) + "_jsaa/mod_1/deg_" + str(deg+1) + "/rep_" + str(data_rep+1) + "/ols/samp_" + str(jack_samp_num[cov_dim][samp_size]+1) + "/jsaa_obj.txt") as wass_file:
					wass_list = wass_file.read().splitlines() 
					for cov_rep in range(num_covariate_replicates):
						JSAA_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] = [float(i) for i in wass_list[cov_rep*num_saa_replicates:(cov_rep+1)*num_saa_replicates]]
						JSAA_gap = JSAA_data[cov_dim,deg,cov_rep,data_rep,samp_size,:] - SAA_data[deg,cov_rep,:]
						JSAA_gap_mean = np.mean(JSAA_gap)
						JSAA_gap_var = np.var(JSAA_gap)
						JSAA_ucb[cov_dim,deg,cov_rep,data_rep,samp_size] = (100.0)*(JSAA_gap_mean + t_value*np.sqrt(JSAA_gap_var/num_saa_replicates))
				wass_file.close()


# Plot ER-SAA + OLS vs J-SAA + OLS and Wasserstein ER-DRO + OLS using Algorithms 2 and 3 for d_x = 100
box_colors = ['limegreen', 'black', 'blue', 'red']
plot_y_max = [70.0, 70.0, 70.0]
for deg in range(num_degrees):
	for cov in range(num_covariate_dimensions):
		plt.figure(figsize=fig_size, dpi=dpi_val)
		plt.rcParams["font.weight"] = "bold"
		plt.rcParams["axes.labelweight"] = "bold"
		plt.rcParams["axes.titleweight"] = "bold"
		plt.rcParams.update({'font.size': font_size})

		for meth in range(4):
			data = np.zeros((num_data_replicates*num_covariate_replicates,num_sample_sizes))
			box_positions = np.zeros(num_sample_sizes)
			box_width = 0.5
			for samp_size in range(num_sample_sizes):
				if meth == 0:
					data[:,samp_size] = JSAA_ucb[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1 + 2.2*samp_size
				elif meth == 1:
					data[:,samp_size] = wass_ERDRO_scrambled_ucb[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 1.5 + 2.2*samp_size
				elif meth == 2:
					data[:,samp_size] = wass_ERDRO_doubleest_ucb[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 2 + 2.2*samp_size
				elif meth == 3:
					data[:,samp_size] = ERSAA_ucb[cov,deg,:,:,samp_size].flatten()
					box_positions[samp_size] = 2.5 + 2.2*samp_size

			plt.boxplot(data,whis=whiskers, showfliers=False, 
						positions=box_positions, widths=box_width,
						boxprops={'linewidth':line_width,'color':box_colors[meth]}, 
						medianprops={'linewidth':line_width,'color':box_colors[meth]},
						whiskerprops={'linewidth':line_width,'color':box_colors[meth]},
						capprops={'linewidth':line_width,'color':box_colors[meth]})
			if deg == 0:
				plt.ylabel('UCB on % optimality gap')
			plt.title(r'$\theta$' + ' = ' + str(degrees[deg]))
			x1,x2,y1,y2 = plt.axis()  
			plt.axis((x1,2.2*(num_sample_sizes)+0.7,0,plot_y_max[deg]))


			box_positions_overall = np.zeros(5*num_sample_sizes)
			box_labels_overall = np.empty(5*num_sample_sizes,dtype=object)
			for samp_size in range(num_sample_sizes):
				box_positions_overall[5*samp_size] = 1 + 2.2*samp_size
				box_labels_overall[5*samp_size] = 'J'
				box_positions_overall[5*samp_size+1] = 1.5 + 2.2*samp_size
				box_labels_overall[5*samp_size+1] = '2'
				box_positions_overall[5*samp_size+2] = 1.75 + 2.2*samp_size
				box_labels_overall[5*samp_size+2] = '\nn=' + str(sample_sizes[cov][samp_size])
				box_positions_overall[5*samp_size+3] = 2 + 2.2*samp_size
				box_labels_overall[5*samp_size+3] = '3'
				box_positions_overall[5*samp_size+4] = 2.5 + 2.2*samp_size
				box_labels_overall[5*samp_size+4] = 'E'
			
			plt.xticks(box_positions_overall, box_labels_overall)
			plt.tick_params(axis='x',
							which='both',
							bottom=False)

			if deg == 0:
				plt.subplots_adjust(left=0.16,right=0.99,bottom=0.12)
			else:
				plt.subplots_adjust(left=0.11,right=0.99,bottom=0.12)
			plt.savefig("figures/fig8_cov" + str(cov+3) + "_deg" + str(deg+1) + ".eps", format='eps')
