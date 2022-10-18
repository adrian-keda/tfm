import pandas as pd
import numpy as np

basedir = '/home/amaqueda/adrian_TFM'

filedir = basedir + '/00_TCGA_annotated_vars_toy_example'

input_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.rare_0.001.'
output_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.rare_0.001.'




def gene_vs_cancer_frequency(df_delmis, df_clinvar, TCGA_study_size):
    
    """
    Function to get the frequencies of individuals who have pathogenic variants per gene and cancer type.

    Parameters:
    df_delmis: dataframe with rare pathogenic variants according to deleterious missense criteria.
    df_clinvar: dataframe with rare pathogenic variants according to clinvar criteria.
    TCGA_study_size:
    """

    # DelMis frequency calculation
    DelMis_unique = df_delmis[['INDIVIDUAL', 'SYMBOL', 'CANCER']].value_counts().reset_index(name = 'count')[['SYMBOL', 'CANCER']].value_counts().reset_index(name = 'count')
    DelMis_count = pd.pivot_table(DelMis_unique, index = 'SYMBOL', columns = 'CANCER', values = 'count', fill_value = 0)
    DelMis_count.insert(0, 'PANCANCER', DelMis_count.sum(axis = 1))

    for i in DelMis_count.columns:
        DelMis_count[i] = DelMis_count[i] / int(TCGA_study_size.loc[i])

    # ClinVar frequency calculation
    ClinVar_unique = df_clinvar[['INDIVIDUAL', 'SYMBOL', 'CANCER']].value_counts().reset_index(name = 'count')[['SYMBOL', 'CANCER']].value_counts().reset_index(name = 'count')
    ClinVar_count = pd.pivot_table(ClinVar_unique, index = 'SYMBOL', columns = 'CANCER', values = 'count', fill_value = 0)
    ClinVar_count.insert(0, 'PANCANCER', ClinVar_count.sum(axis = 1))

    for i in ClinVar_count.columns:
        ClinVar_count[i] = ClinVar_count[i] / int(TCGA_study_size.loc[i])
    
    return DelMis_count, ClinVar_count




# Load dataframes
DelMis_variants = pd.read_csv(input_file + 'DelMis.chr21.tsv', sep = '\t', header = 0)
ClinVar_variants = pd.read_csv(input_file + 'ClinVar.chr21.tsv', sep = '\t', header = 0)
TCGA_study_size = pd.read_csv(basedir + '/TCGA_cancer_count.tsv', sep = '\t', header = 0, index_col = 0)


# Calculate frequencies
DelMis_freq, ClinVar_freq = gene_vs_cancer_frequency(DelMis_variants, ClinVar_variants, TCGA_study_size)


# Save dataframes
DelMis_freq.to_csv(path_or_buf = output_file + 'DelMis.frequencies.chr21.tsv', sep = '\t', index = True)
ClinVar_freq.to_csv(path_or_buf = output_file + 'ClinVar.frequencies.chr21.tsv', sep = '\t', index = True)