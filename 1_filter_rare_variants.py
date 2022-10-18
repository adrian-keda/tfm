import pandas as pd
import numpy as np

basedir = '/home/amaqueda/adrian_TFM'

filedir = basedir + '/00_TCGA_annotated_vars_toy_example'

thr = 0.001

input_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.chr21.tsv'
output_file = filedir + f'/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.rare_{thr}.chr21.tsv'




def rare_variant_filter(df, threshold = 0.01):

    """
    This functions filters rare variants from a dataframe containing frequency information of variants on TCGA and gnomAD.
    Variants whose frequencies are lower than the specified threshold are retained. By default, the threshold is 0.01 (rare variants).

    Parameters:
        df: dataframe containing information about variants.
        threshold: user specified threshold.
    """

    a = df[pd.isnull(df['gnomAD_AF']) & (df['TCGA_ALL_AF'] <= threshold) & (df['TCGA_EUR_AF'] <= threshold) & (df['TCGA_AMR_AF'] <= threshold)
    & (df['TCGA_EAS_AF'] <= threshold) & (df['TCGA_AFR_AF'] <= threshold)]

    b = df[(df['gnomAD_AF'] <= threshold) & (df['gnomAD_AFR_AF'] <= threshold) & (df['gnomAD_AMR_AF'] <= threshold) & (df['gnomAD_EAS_AF'] <= threshold)
    & (df['gnomAD_SAS_AF'] <= threshold) & (df['gnomAD_NFE_AF'] <= threshold) & (df['TCGA_ALL_AF'] <= threshold) & (df['TCGA_EUR_AF'] <= threshold)
    & (df['TCGA_AMR_AF'] <= threshold) & (df['TCGA_EAS_AF'] <= threshold) & (df['TCGA_AFR_AF'] <= threshold)]

    rare_variants = pd.concat([a, b]).reset_index(drop = True).drop_duplicates()

    return rare_variants




# Load dataframe
df = pd.read_csv(input_file, sep = '\t', header = 0, na_values = '-')


# Get variants whose frequency for every population is under the specified threshold.
df = rare_variant_filter(df, threshold = thr)


# Save dataframe
df.to_csv(path_or_buf = output_file, sep = '\t', index = False)