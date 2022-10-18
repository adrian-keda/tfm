import pandas as pd
import numpy as np

basedir = '/home/amaqueda/adrian_TFM'




def add_pop_cancer (df, pop_cancer = pd.read_csv(basedir + '/barcode_study_ancestry.tsv', sep = '\t', header = 0)):

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




def calc_frequencies_pop_2(df, individual_count = pd.read_csv(basedir + '/TCGA_ancestry_count.tsv', sep = '\t', header = 0), pops = ['ALL', 'EUR', 'AMR', 'EAS', 'AFR']):

    """
    Function used to calculate the frequencies of each variant in the desired populations.

    Parameters:
    df: dataframe with information of variants.
    individual_count: dataframe containing the number of individuals for each population.
    pops: list to specify in which populations calculate the frequencies.
    """

    variants = df['Otherinfo3'].unique()   # list of variants
    FREQS = np.zeros(shape = (len(variants), len(pops)) )  # array to store frequencies

    for pop in pops:
        row = 0    # fila del array
        if pop == 'ALL':
            df_pop = df
            n_individuals = int(individual_count[individual_count['ancestry'] == 'TOTAL']['count']) # Number of individuals from all populations
        else:    
            df_pop = df[df['POPULATION'] == pop]
            n_individuals = int(individual_count[individual_count['ancestry'] == pop]['count']) # Number of individuals from a single population
        
        col = pops.index(pop)   # columna del array
        a = df_pop[['Otherinfo3', 'GT']].value_counts().reset_index(name='count')
        for var in variants:

            # HETEROCIGOTOS
            try:
                hetz = int(a[ (a['Otherinfo3'] == var) & (a['GT'] == '0/1') ]['count'])
            except TypeError:
                hetz = 0
            
            # HOMOCIGOTOS
            try:
                homoz = int(a[ (a['Otherinfo3'] == var) & (a['GT'] == '1/1') ]['count']) * 2
            except TypeError:
                homoz = 0
            
            FREQS[row, col] = (homoz + hetz) / (n_individuals * 2)  # Computing allele frequency
            row = row + 1
    
    data = {'Otherinfo3':variants.tolist(), 'TCGA_ALL_AF':FREQS[:,0], 'TCGA_EUR_AF':FREQS[:,1], 
            'TCGA_AMR_AF':FREQS[:,2], 'TCGA_EAS_AF':FREQS[:,3], 'TCGA_AFR_AF':FREQS[:,4]}
    df_frequencies = pd.DataFrame(data = data)

    return df_frequencies




def add_allele_population_frequencies(df, pop_allele_frequencies):
    
    """
    Function used to add the frequencies of each variant on each population of the study.

    Parameters:
    df: dataframe with information of variants.
    pop_allele_frequencies: dataframe with frequencies of variants on different populations.
    """
    
    vars = list(df['Otherinfo3'])
    vars_unique = list(pop_allele_frequencies['Otherinfo3'])
    all_af = list(pop_allele_frequencies['TCGA_ALL_AF'])
    eur_af = list(pop_allele_frequencies['TCGA_EUR_AF'])
    amr_af = list(pop_allele_frequencies['TCGA_AMR_AF'])
    eas_af = list(pop_allele_frequencies['TCGA_EAS_AF'])
    afr_af = list(pop_allele_frequencies['TCGA_AFR_AF'])

    frec_dict = {}
    for var in vars_unique:
        pos = vars_unique.index(var)
        frec_dict[var] = {'ALL':all_af[pos], 'EUR':eur_af[pos], 'AMR':amr_af[pos], 'EAS':eas_af[pos], 'AFR':afr_af[pos]}
    
    ALL_AF = []
    EUR_AF = []
    AMR_AF = []
    EAS_AF = []
    AFR_AF = []

    for var in vars:
        ALL_AF.append(frec_dict[var]['ALL'])
        EUR_AF.append(frec_dict[var]['EUR'])
        AMR_AF.append(frec_dict[var]['AMR'])
        EAS_AF.append(frec_dict[var]['EAS'])
        AFR_AF.append(frec_dict[var]['AFR'])
    
    df['TCGA_ALL_AF'] = ALL_AF
    df['TCGA_EUR_AF'] = EUR_AF
    df['TCGA_AMR_AF'] = AMR_AF
    df['TCGA_EAS_AF'] = EAS_AF
    df['TCGA_AFR_AF'] = AFR_AF

    return df




filedir = basedir + '/00_TCGA_annotated_vars_toy_example'

input_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.chr21.tsv'
output_file = filedir + '/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.chr21.tsv'


# Read the .tsv file
df = pd.read_csv(input_file, sep='\t', header = 0)


# Add INDIVIDUAL and SAMPLE_NUMERIC_CODE columns
df['INDIVIDUAL'] = [i[0:12] for i in df['SAMPLE']]
df['SAMPLE_NUM_CODE'] = [i[13:15] for i in df['SAMPLE']]


# Remove fully duplicated rows
df = df.drop_duplicates()


# Remove rows dupplicated due to different annotation
df = df.drop_duplicates(['Otherinfo3', 'INDIVIDUAL', 'SAMPLE_NUM_CODE'], keep='first')


# Add POPULATION and CANCER_TYPE columns
df = add_pop_cancer(df)


# Calculating variants' frequencies on five populations (ALL, EUR, AMR, EAS, AFR)
df_frequencies = calc_frequencies_pop_2(df[['Otherinfo3', 'GT', 'INDIVIDUAL', 'POPULATION']])


# Adding allele frequencies
df = add_allele_population_frequencies(df, df_frequencies)

# Saving dataframe for future analysis
df.to_csv(path_or_buf = output_file, sep = '\t', index = False)