import maps_processing 
import confidence_calculation 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Path to files
file_names = ['test_data/test_data_rep1.cool','test_data/test_data_rep2.cool']
tads_positions = pd.read_csv('test_data/tads.csv')
tads_positions.columns = ['Chr','Start','End']

# Preliminary step to choose the best parameter for normalization
k_parameter = maps_processing.choose_k(file_names)

# Generate confidence for two Hi-C matrices
mtx1, regions, bin_size = confidence_calculation.download_matrix(file_names[0])
mtx2, regions, bin_size = confidence_calculation.download_matrix(file_names[1])
chrm = regions[0] # Choose chromosome
confidence = confidence_calculation.calculate_confidence(mtx1(chrm,chrm),mtx2(chrm,chrm),k_parameter)

# Generate average normalized values for areas around tads 
avg_tad = maps_processing.average_tad(file_names, k_parameter, tads_positions, area_size=7)
sns.heatmap(avg_tad['rep_1'],cmap='RdBu_r', vmin=0.8, square=True, xticklabels=False, yticklabels=False)
plt.show()

# Generate normalized maps by downsampling
maps_processing.downsampling(file_names, k_parameter, factor=2, output_file_name='output')



