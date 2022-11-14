import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]




def add_pop_cancer (df, pop_cancer = pd.read_csv(snakemake.input[1], sep = '\t', header = 0)):

    """
    Function used to add the populations and cancer type columns to the dataframe.

    Parameters:
    df: dataframe with information of variants.
    pop_cancer: dataframe
    """

    ind_list = list(pop_cancer['TCGA.Barcode'])
    pop_list = list(pop_cancer['TCGA.Ancestry'])
    cancer_list = list(pop_cancer['TCGA.Study.Cancer'])

    dict = {}
    for ind in ind_list:
        dict[ind] = [pop_list[ind_list.index(ind)], cancer_list[ind_list.index(ind)]]

    ind_list_var = list(df['INDIVIDUAL'])
    pop_list_var = []
    cancer_list_var = []
    
    for ind in ind_list_var:
        pop_list_var.append(dict[ind][0])
        cancer_list_var.append(dict[ind][1])

    df['POPULATION'] = pop_list_var
    df['CANCER'] = cancer_list_var

    return df




def calc_frequencies_pop(df, individual_count = pd.read_csv('/home/amaqueda/adrian_TFM/TCGA_ancestry_count.tsv', sep = '\t', header = 0)):

    """
    Function used to calculate the frequencies of each variant in the desired populations.

    Parameters:
    df: dataframe with information of variants.
    individual_count: dataframe containing the number of individuals for each population.
    """
    
    variants = df['Otherinfo3'].unique() # list of variants
    FREQS = []  # list to store frequencies
    n_individuals = individual_count['count'].sum() # number of individuals in or samples

    a = df[['Otherinfo3', 'GT']].value_counts().reset_index(name='count')
    a['GT'] = [1 if i == '0/1' else 2 for i in a['GT']]

    for var in variants:
        
        b = a[ a['Otherinfo3'] == var ]
        n_alleles = (b['count'] * b['GT']).sum() # Compute number of alleles
        
        FREQS.append((n_alleles) / (n_individuals * 2)) # Computing and storing allele frequency
        
    return variants, FREQS




# Read the tsv file
df = pd.read_csv(input_file, sep='\t', header = 0)


# Add INDIVIDUAL and SAMPLE_NUMERIC_CODE columns
df['INDIVIDUAL'] = [i[0:12] for i in df['SAMPLE']]
df['SAMPLE_NUM_CODE'] = [i[13:15] for i in df['SAMPLE']]



# Remove rows dupplicated due to different annotation
df = df.drop_duplicates(['Otherinfo3', 'INDIVIDUAL', 'SAMPLE_NUM_CODE'], keep='first')


# Add POPULATION and CANCER_TYPE columns
df = add_pop_cancer(df)


# Calculating variants' frequencies
variants, frequencies = calc_frequencies_pop(df[['Otherinfo3', 'GT', 'INDIVIDUAL', 'POPULATION']])
variants = variants.tolist()


# Creating a dictionary variant : frequency
dict = {}
for var in variants:
    dict[var] = frequencies[variants.index(var)]

# Adding allele frequencies
df['TCGA_ALL_AF'] = [ dict[variant] for variant in df['Otherinfo3'] ]


# Saving dataframe for future analysis
df.to_csv(path_or_buf = output_file, sep = '\t', index = False)