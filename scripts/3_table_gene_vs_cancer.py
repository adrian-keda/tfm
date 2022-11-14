import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]




def gene_vs_cancer_frequency_DelMis(df, TCGA_study_size):
    
    """
    Function to get the frequencies of rare pathogenic variants according to DelMis criteria per gene and cancer type.

    Parameters:
    df: dataframe with rare pathogenic variants.
    TCGA_study_size: dataframe with the number of individuals per cancer type, including pancancer.
    """

    DelMis_unique = df[['INDIVIDUAL', 'SYMBOL', 'CANCER']].value_counts().reset_index(name = 'vars_per_ind')[['SYMBOL', 'CANCER']].value_counts().reset_index(name = 'n_inds')
    DelMis_count = pd.pivot_table(DelMis_unique, index = 'SYMBOL', columns = 'CANCER', values = 'n_inds', fill_value = 0)

    for i in TCGA_study_size.index:
        if i not in DelMis_count.columns:
            DelMis_count[i] = [0] * DelMis_count.shape[0]


    DelMis_count['PANCANCER'] = DelMis_count.sum(axis = 1)

    for i in DelMis_count.columns:
        DelMis_count[i] = DelMis_count[i] / int(TCGA_study_size.loc[i])
    
    return DelMis_count




def gene_vs_cancer_frequency_ClinVar(df, TCGA_study_size):
    
    """
    Function to get the frequencies of rare pathogenic variants according to ClinVar criteria per gene and cancer type.

    Parameters:
    df: dataframe with rare pathogenic variants.
    TCGA_study_size: dataframe with the number of individuals per cancer type, including pancancer.
    """

    ClinVar_unique = df[['INDIVIDUAL', 'SYMBOL', 'CANCER']].value_counts().reset_index(name = 'vars_per_ind')[['SYMBOL', 'CANCER']].value_counts().reset_index(name = 'n_inds')
    ClinVar_count = pd.pivot_table(ClinVar_unique, index = 'SYMBOL', columns = 'CANCER', values = 'n_inds', fill_value = 0)

    for i in TCGA_study_size.index:
        if i not in ClinVar_count.columns:
            ClinVar_count[i] = [0] * ClinVar_count.shape[0]


    ClinVar_count['PANCANCER'] = ClinVar_count.sum(axis = 1)

    for i in ClinVar_count.columns:
        ClinVar_count[i] = ClinVar_count[i] / int(TCGA_study_size.loc[i])
    
    return ClinVar_count




# Load dataframes
df = pd.read_csv(input_file, sep = '\t', header = 0)
TCGA_study_size = pd.read_csv(snakemake.input[1], sep = '\t', header = 0, index_col = 0)

# Calculate DelMis variants frequencies in cancer
if 'DelMis' in input_file:
    DelMis_freq = gene_vs_cancer_frequency_DelMis(df, TCGA_study_size)
    DelMis_freq.to_csv(path_or_buf = output_file, sep = '\t', index = True)


# Calculate ClinVar variants frequencies in cancer
if 'ClinVar' in input_file:
    ClinVar_freq = gene_vs_cancer_frequency_ClinVar(df, TCGA_study_size)
    ClinVar_freq.to_csv(path_or_buf = output_file, sep = '\t', index = True)