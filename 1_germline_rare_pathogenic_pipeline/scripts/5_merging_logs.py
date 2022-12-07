import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




# Creating lists of dataframes for chromosomes' filters logs and pathogenic filters logs
chr_logs_dfs = [pd.read_csv(file, sep = '\t', header = 0) for file in snakemake.input["chr_logs"]]
patho_logs_dfs = [pd.read_csv(file, sep = '\t', header = 0) for file in snakemake.input["patho_logs"]]


# Creating dataframe combining the information for each field in the chromosomes' filters logs. It's just a sum of fields  :)
combined_steps_df = pd.DataFrame(data = {
    'Filter_ID' : chr_logs_dfs[0]['Filter_ID'].tolist(),
    'Input' : np.array([log['Input'] for log in chr_logs_dfs]).sum(axis = 0),
    'Output' : np.array([log['Output'] for log in chr_logs_dfs]).sum(axis = 0)
})


# Merging dataframe with the number of variants before and after each filter for all the chromosomes and the
# dataframes for the pathogenic filters.
final_df = pd.concat([combined_steps_df] + [i for i in patho_logs_dfs], axis = 0, ignore_index = False)
final_df['Difference'] = final_df['Input'] - final_df['Output']

# Saving the dataframe
final_df.to_csv(path_or_buf = snakemake.output["output_log"], sep = '\t', index = False)