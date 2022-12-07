import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input["input"]
output_file = snakemake.output["output"]
vars_filtered_out = snakemake.output["vars_filtered_log"]



def filter_pathogenic_variants_DelMis(df):

    """
    Function to filter rare pathogenic variants according to DelMis criteria (deleterious missense variants). We use MetaLR predictor.

    Parameters:
    df: dataframe with information of variants.
    """

    #Filter missense variants first, then deleterious missense ones.
    mis_variants = df[(df['Consequence'] == 'missense_variant') | (df['Consequence'] == 'missense_variant,splice_region_variant') | (df['Consequence'] == 'missense_variant,NMD_transcript_variant')]

    del_mis_variants = mis_variants[mis_variants['MetaLR_pred'] == 'D']
    
    return del_mis_variants




def filter_pathogenic_variants_ClinVar(df,
            patho_list_in = ['risk_factor','pathogenic','likely_pathogenic','likely_pathogenic_risk_factor', 'pathogenic/likely_pathogenic', '_risk_factor'],
            patho_list_out = ['conflicting_interpretations_of_pathogenicity', 'uncertain_significance','not_provided', 'benign', 'likely_benign']):

    """
    Function to filter rare pathogenic variants according to ClinVar criteria (clinvar significance variants).

    Parameters:
    df: dataframe with information of variants.
    path_list: criteria to filter pathogenic variants.
    """

    # Filter variants which have at least one of the ClinSig tags in "patho_list_in"
    picklist = []
    for index, row in df.iterrows():
        picklist.append(any(x in patho_list_in for x in str(row['CLIN_SIG']).split(',')))
    df = df.loc[picklist]

    # Filter out variants which have at least one of the ClinSig tags in "patho_list_out"
    picklist = []
    for index, row in df.iterrows():
        picklist.append(not any(x in patho_list_out for x in str(row['CLIN_SIG']).split(','))) 
    clin_variants = df.loc[picklist]
    
    return clin_variants




# Load dataframe
df = pd.read_csv(input_file, sep = '\t', header = 0)
max_vars = df.shape[0]


# Filter pathogenic variants according to DelMis criteria
if 'DelMis' in output_file:
    DelMis_variants = filter_pathogenic_variants_DelMis(df)
    DelMis_variants.to_csv(path_or_buf = output_file, sep = '\t', index = False)
    
    DelMis_vars = DelMis_variants.shape[0]
    new_row = {'Filter_ID':['DelMis'],
            'Input':[max_vars],
            'Output':[DelMis_vars]}
    DelMis_log = pd.DataFrame(data = new_row)
    DelMis_log.to_csv(path_or_buf = vars_filtered_out, sep = '\t', index = False)


# Filter pathogenic variants according to ClinVar criteria
if 'ClinVar' in output_file:
    ClinVar_variants = filter_pathogenic_variants_ClinVar(df)
    ClinVar_variants.to_csv(path_or_buf = output_file, sep = '\t', index = False)
    
    ClinVar_vars = ClinVar_variants.shape[0]
    new_row = {'Filter_ID':['ClinVar'],
            'Input':[max_vars],
            'Output':[ClinVar_vars]}
    ClinVar_log = pd.DataFrame(data = new_row)
    ClinVar_log.to_csv(path_or_buf = vars_filtered_out, sep = '\t', index = False)