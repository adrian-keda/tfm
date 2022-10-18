import pandas as pd
import numpy as np

basedir = '/home/amaqueda/adrian_TFM'

filedir = basedir + '/00_TCGA_annotated_vars_toy_example'

input_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.rare_0.001.chr21.tsv'
output_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.rare_0.001.'




def filter_pathogenic_variants(df):

    """
    Function to filter rare pathogenic variants according to two criteria (deleterious missense variants and clinvar).
    To select deleterious variants we use MetaLR predictor.

    Parameters:
    df: dataframe with information of variants.
    """

    df[['Ref', 'Alt']] = df[['Ref', 'Alt']].fillna('-')

    #Filter missense variants first, then deleterious missense ones.
    mis_variants = df[(df['Consequence'] == 'missense_variant') | (df['Consequence'] == 'missense_variant,splice_region_variant')]

    del_mis_variants = mis_variants[mis_variants['MetaLR_pred'] == 'D']

    #Filter pathogenic variants acordint to ClinVar.
    clin_variants = df[(df['CLIN_SIG'] == 'pathogenic') | (df['CLIN_SIG'] == 'pathogenic/likely_pathogenic') | 
                   (df['CLIN_SIG'] == 'risk_factor') | (df['CLIN_SIG'] == 'likely_pathogenic_risk_factor')]
    
    return del_mis_variants, clin_variants




# Load dataframe
df = pd.read_csv(input_file, sep = '\t', header = 0)


# Filter pathogenic variants according to ClinVar and DelMis criteria
DelMis_variants, ClinVar_variants = filter_pathogenic_variants(df)


# Save dataframes
DelMis_variants.to_csv(path_or_buf = output_file + 'DelMis.chr21.tsv', sep = '\t', index = False)
ClinVar_variants.to_csv(path_or_buf = output_file + 'ClinVar.chr21.tsv', sep = '\t', index = False)