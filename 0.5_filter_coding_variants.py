import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]




def filter_coding_variants (df):

    out_df = df[(df['Func.refGene'] == 'exonic') | (df['Func.refGene'] == 'exonic;splicing') | (df['Func.refGene'] == 'splicing') |
                (df['Func.refGene'] == 'UTR3') | (df['Func.refGene'] == 'UTR5') | (df['Func.refGene'] == 'UTR5;UTR3')]

    return out_df 




# Load dataframe
df = pd.read_csv(input_file, sep = '\t', header = 0)


# Filter variants which affect the canonical isotope
df = df[df['CANONICAL'] == 'YES']


# Filter variants by their consequence
list_consequence_interest = ['missense_variant', 'synonymous_variant', '3_prime_UTR_variant', '5_prime_UTR_variant',
       'splice_region_variant,synonymous_variant', 'splice_region_variant,intron_variant', 'frameshift_variant',
       'missense_variant,splice_region_variant', 'inframe_deletion','splice_acceptor_variant', 'stop_gained,splice_region_variant',
       'stop_gained', 'frameshift_variant,splice_region_variant', 'splice_donor_variant',
       'splice_donor_variant,non_coding_transcript_variant','start_lost', 'splice_region_variant,non_coding_transcript_exon_variant',
       'splice_acceptor_variant,non_coding_transcript_variant', 'stop_retained_variant', 'splice_acceptor_variant,intron_variant',
       'splice_region_variant,3_prime_UTR_variant', 'inframe_deletion,splice_region_variant',
       'splice_donor_variant,coding_sequence_variant', 'stop_gained,frameshift_variant', 'stop_lost',
       'splice_region_variant,intron_variant,non_coding_transcript_variant', 'splice_donor_variant,intron_variant',
       'splice_donor_variant,coding_sequence_variant,intron_variant', 'splice_region_variant,5_prime_UTR_variant',
       'stop_gained,inframe_insertion', 'inframe_insertion,splice_region_variant',
       'stop_gained,inframe_deletion', 'start_lost,splice_region_variant', 'frameshift_variant,stop_lost',
       'frameshift_variant,start_lost,start_retained_variant']

consequence = df['Consequence'].tolist()
bool_list = [True if i in list_consequence_interest else False for i in consequence] # Boolean list to pick rows later
df = df.loc[bool_list]

# Filter variants in coding regions
out_df = filter_coding_variants (df)


# Save
out_df.to_csv(path_or_buf = output_file, sep = '\t', index = False)