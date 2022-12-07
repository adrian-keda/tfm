import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input["input"]
vars_filtered_log_in = snakemake.input["vars_filtered_log"]

output_file = snakemake.output["output"]
vars_filtered_log_out = snakemake.output["vars_filtered_log"]



def filter_coding_variants (df):

    out_df = df[(df['Func.refGene'] == 'exonic') | (df['Func.refGene'] == 'exonic;splicing') | (df['Func.refGene'] == 'splicing') |
                (df['Func.refGene'] == 'UTR3') | (df['Func.refGene'] == 'UTR5') | (df['Func.refGene'] == 'UTR5;UTR3')]

    return out_df 




# Load dataframes
df = pd.read_csv(input_file, sep = '\t', header = 0)
variant_log = pd.read_csv(vars_filtered_log_in, sep = '\t', header = 0)


# Filter variants which affect the canonical isotope
max_vars = df.shape[0]
df = df[df['CANONICAL'] == 'YES']
canonical_vars = df.shape[0]


# Filter variants by their consequence
list_consequence_interest = ['missense_variant', 'synonymous_variant', '3_prime_UTR_variant', '5_prime_UTR_variant',
       'splice_region_variant,synonymous_variant', 'frameshift_variant',
       'missense_variant,splice_region_variant', 'inframe_deletion','splice_acceptor_variant', 'stop_gained,splice_region_variant',
       'stop_gained', 'frameshift_variant,splice_region_variant', 'splice_donor_variant', 'start_lost',
       'stop_retained_variant', 'splice_region_variant,3_prime_UTR_variant',
       'inframe_deletion,splice_region_variant', 'splice_donor_variant,coding_sequence_variant', 'stop_gained,frameshift_variant', 'stop_lost',
       'splice_region_variant,5_prime_UTR_variant', 'stop_gained,inframe_insertion',
       'inframe_insertion,splice_region_variant', 'stop_gained,inframe_deletion', 'start_lost,splice_region_variant',
       'frameshift_variant,stop_lost', 'frameshift_variant,start_lost,start_retained_variant']

consequence = df['Consequence'].tolist()
bool_list = [True if i in list_consequence_interest else False for i in consequence] # Boolean list to pick rows later
df = df.loc[bool_list]
consequence_vars = df.shape[0]


# Filter variants in coding regions
out_df = filter_coding_variants (df)
coding_vars = out_df.shape[0]


# Data for log file
new_rows = {'Filter_ID':['canonical', 'consequence', 'coding'],
            'Input':[max_vars, canonical_vars, consequence_vars],
            'Output':[canonical_vars, consequence_vars, coding_vars]}
variant_log_new_rows = pd.DataFrame(data = new_rows)
variant_log = pd.concat([variant_log, variant_log_new_rows], axis = 0, ignore_index = False)

# Save
out_df.to_csv(path_or_buf = output_file, sep = '\t', index = False)
variant_log.to_csv(path_or_buf = vars_filtered_log_out, sep = '\t', index = False)