import pandas as pd
import numpy as np




# Defining parameters 
thr = snakemake.params[0]
input_file = snakemake.input["input"]
vars_filtered_log_in = snakemake.input["vars_filtered_log"]

output_file = snakemake.output["output"]
vars_filtered_log_out = snakemake.output["vars_filtered_log"]




# def rare_variant_filter(df, threshold = 0.01):

#     """
#     This functions filters rare variants from a dataframe containing frequency information of variants on TCGA and gnomAD.
#     Variants whose frequencies are lower than the specified threshold are retained. By default, the threshold is 0.01 (rare variants).

#     Parameters:
#         df: dataframe containing information about variants.
#         threshold: user specified threshold.
#     """
#     #Changing - to NaN in gnomAD's allele frequencies columns
#     df[['gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_SAS_AF','gnomAD_NFE_AF']]=df[['gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_SAS_AF','gnomAD_NFE_AF']].replace('-',np.nan)
    
#     # Filtering step for variants with no allele frequencies in gnomAD
#     a = df[pd.isnull(df['gnomAD_AF']) & (df['TCGA_ALL_AF'] <= threshold)]
    
#     # Creating dataframe for variants with allele frequencies in gnomAD
#     b = df[[not i for i in pd.isnull(df['gnomAD_AF'])]]
    
#     # Changing column type from str to float
#     b = b.astype({'gnomAD_AF' : 'float', 'gnomAD_AFR_AF' : 'float', 'gnomAD_AMR_AF': 'float', 'gnomAD_EAS_AF' : 'float', 'gnomAD_SAS_AF' : 'float', 'gnomAD_NFE_AF' : 'float'})
    
#     # Filtering step for variants with allele frequencies in gnomAD
#     b = b[(b['gnomAD_AF'] <= threshold) & (b['gnomAD_AFR_AF'] <= threshold) & (b['gnomAD_AMR_AF'] <= threshold) & (b['gnomAD_EAS_AF'] <= threshold) & (b['gnomAD_SAS_AF'] <= threshold) & (b['gnomAD_NFE_AF'] <= threshold) & (b['TCGA_ALL_AF'] <= threshold)]
    
#     # Concat dataframes
#     rare_variants = pd.concat([a, b]).reset_index(drop = True).drop_duplicates()

#     return rare_variants




def rare_variant_filter(df, threshold = 0.01):

    """
    This functions filters rare variants from a dataframe containing frequency information of variants on TCGA and gnomAD.
    Variants whose frequencies are lower than the specified threshold are retained. By default, the threshold is 0.01 (rare variants).

    Parameters:
        df: dataframe containing information about variants.
        threshold: user specified threshold.
    """

    df[['gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_SAS_AF','gnomAD_NFE_AF']]=df[['gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_EAS_AF','gnomAD_SAS_AF','gnomAD_NFE_AF']].replace('-',np.nan)

    a = df[pd.isnull(df['gnomAD_AF']) & (df['TCGA_ALL_AF'] <= threshold)]
    
    b = df.loc[[not i for i in pd.isnull(df['gnomAD_AF'])]]
    
    b = b.astype({'gnomAD_AF' : 'float', 'gnomAD_AFR_AF' : 'float', 'gnomAD_AMR_AF': 'float', 'gnomAD_EAS_AF' : 'float', 'gnomAD_SAS_AF' : 'float', 'gnomAD_NFE_AF' : 'float'})

    b = b[(b['gnomAD_AF'] <= threshold) & (b['gnomAD_AFR_AF'] <= threshold) & (b['gnomAD_AMR_AF'] <= threshold) & (b['gnomAD_EAS_AF'] <= threshold) & (b['gnomAD_SAS_AF'] <= threshold) & (b['gnomAD_NFE_AF'] <= threshold) & (b['TCGA_ALL_AF'] <= threshold)]

    rare_variants = pd.concat([a, b]).reset_index(drop = True).drop_duplicates()

    return rare_variants




# Load dataframe
df = pd.read_csv(input_file, sep = '\t', header = 0, na_values = '-')
variant_log = pd.read_csv(vars_filtered_log_in, sep = '\t', header = 0)
max_vars = df.shape[0]


# Get variants whose frequency for every population is under the specified threshold.
df = rare_variant_filter(df, threshold = thr)
rare_vars = df.shape[0]


# Data for log file
new_row = {'Filter_ID':['rare'],
            'Input':[max_vars],
            'Output':[rare_vars]}
new_row = pd.DataFrame(data = new_row)
variant_log = pd.concat([variant_log, new_row], axis = 0, ignore_index = True)


# Save dataframes
df.to_csv(path_or_buf = output_file, sep = '\t', index = False)
variant_log.to_csv(path_or_buf = vars_filtered_log_out, sep = '\t', index = False)