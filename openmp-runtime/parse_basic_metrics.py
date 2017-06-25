import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

exp_name = {
'1' : 'blackscholes-omp_p',
'2' : 'bodytrack-omp_p',
'3' : 'freqmine_p',
'4' : 'bt_npb',
'7' : 'ep_npb',
'12' : 'semphy_mb',
'13' : 'rnaseq'
}
exp_filename = {}
exp_df = {}
exp_sf_ranges = {}


for file in os.listdir("."):
	if file.endswith(".metrics") and 'big' in file:
		n_exp = file.split('_')[2]
		if n_exp in exp_name.keys():
			exp = exp_name[n_exp]
			exp_filename[exp] = file

			df_big = pd.read_csv(file)
			df_little = pd.read_csv(file.replace("big","little"))

			del df_big['nsample']
			del df_big['pid']
			del df_big['event']
			del df_little['nsample']
			del df_little['pid']
			del df_little['event']

			df_big.columns = ["B_mips","B_btr", "B_llc_miss_pki"]
			df_little.columns = ["S_mips","S_btr", "S_llc_miss_pki"]


			if len(df_little) == len(df_big):
				merged_df = df_big.copy()
				merged_df['S_mips'] = df_little.S_mips
				merged_df['S_btr'] = df_little.S_btr
				merged_df['S_llc_miss_pki'] = df_little.S_llc_miss_pki
				merged_df['SF'] = merged_df.B_mips / merged_df.S_mips

				exp_df[exp] = merged_df.copy()
				if len(df_little) > 1:
					merged_df.SF.plot()
					plt.savefig("sf_over_time_%s.pdf"%exp)
					plt.close()

				sf_max = merged_df.SF.max()
				sf_min = merged_df.SF.min()
				dif = sf_max - sf_min
				nsteps = np.ceil(dif/0.5)
				total_samples = len(merged_df)


				# Calculate range to have 0.5 partitions starting just before SF_min
				sf_i = sf_min - (sf_min % 0.5)
				sf_ranges = []
				sf_range_names = []
				print "\n%s: min(%f) max(%f) total_samples(%d)" % (exp, sf_min, sf_max, total_samples)

				for i in np.arange(nsteps):
					step_start = sf_i
					step_end = sf_i + 0.5

					mask = (step_start <= merged_df['SF']) & (merged_df['SF'] < step_end)
					sf_samples = merged_df[mask]
					sf_partition_count = len(sf_samples)
					sf_partition_perc = (float(sf_partition_count) / total_samples) * 100

					print "SF_i(%f): %f per100, count(%f)" %(sf_i, sf_partition_perc, sf_partition_count)

					sf_ranges.append(sf_partition_perc)
					sf_range_names.append( "%s-%s"%(str(step_start),str(step_end)) )
					sf_i += 0.5

				ranges_df = pd.DataFrame([],columns=['range','%'])
				ranges_df['range'] = sf_range_names
				ranges_df['%'] = sf_ranges

				ranges_df.to_csv("sf_partitions_%s.csv"%exp, index=False)
				if len(df_little) > 1:
					ax = ranges_df.plot(x='range',y='%',kind='bar')
					plt.rcParams['xtick.major.pad'] = 5
					plt.axis("tight")
					plt.xlabel("Speedup range")
					ax.set_xticklabels(ranges_df.ix[:,0].values,rotation=-45,ha='center')
					plt.savefig("sf_partitions_%s.pdf"%exp)
					plt.close()