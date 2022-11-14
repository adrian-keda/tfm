import pandas as pd
import numpy as np




# Defining parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]




def CPGs_subset (freqs, CPG_metadata):
    
    """
    Function to extract the subset of frequencies from CPGs and nonCPGs

    Parameters:
    -freqs: table with frequencies on genes vs cancer type.
    -CPG_metadata: table with list of CPGs and CPG-like genes.
    """

    CPGs = CPG_metadata[ CPG_metadata ['CPG'] == 'Yes' ]['Gene'].tolist()
    freqs_index = freqs.index.tolist()

    # Loop to get which CPGs from the metadata are present in our frequency table
    CPGs_in_index = []
    for i in CPGs:
        if i in freqs_index:
            CPGs_in_index.append(i)

    # Get subset
    CPGs_subset = freqs.loc[CPGs_in_index]
    
    return CPGs_subset




def nonCPGs_subset (freqs, CPG_metadata):
    
    """
    Function to extract the subset of frequencies from CPGs and nonCPGs

    Parameters:
    -freqs: table with frequencies on genes vs cancer type.
    -CPG_metadata: table with list of CPGs and CPG-like genes.
    """

    CPGs = CPG_metadata[ CPG_metadata ['CPG'] == 'Yes' ]['Gene'].tolist()
    freqs_index = freqs.index.tolist()

    # Loop to get which CPGs from the metadata are present in our frequency table
    CPGs_in_index = []
    for i in CPGs:
        if i in freqs_index:
            CPGs_in_index.append(i)

    # Get genes that are CPGs or CPG-like
    CPGs_and_CPG_like = CPG_metadata['Gene'].tolist()
    
    # Loop to get nonCPGs in index
    non_CPGs_in_index = []
    for i in freqs_index:
        if i not in CPGs_and_CPG_like:
            non_CPGs_in_index.append(i)

    # Get subset
    non_CPGs_subset = freqs.loc[non_CPGs_in_index].sample(len(CPGs_in_index), random_state = 14)
    
    return non_CPGs_subset




## THIS IS JUST IN CASE I MAKE IT TO PRODUCE BOTH FILES AT THE SAME TIME

# def CPGs_nonCPGs_subset (freqs, CPG_metadata):
    
#     CPGs = CPG_metadata[ CPG_metadata ['CPG'] == 'Yes' ]['Gene'].tolist()
#     freqs_index = freqs.index.tolist()

#     # Loop to get which CPGs from the metadata are present in our frequency table
#     CPGs_in_index = []
#     for i in CPGs:
#         if i in freqs_index:
#             CPGs_in_index.append(i)

#     # Get subset
#     CPGs_subset = freqs.loc[CPGs_in_index]


#     # Loop to get which genes from the table are not CPGs or CPG-like genes (basically,they are nonCPGs)
#     CPGs_and_CPG_like = CPG_metadata['Gene'].tolist()
#     non_CPGs_in_index = []

#     for i in freqs_index:
#         if i not in CPGs_and_CPG_like:
#             non_CPGs_in_index.append(i)

#     # Get subset
#     non_CPGs_subset = freqs.loc[non_CPGs_in_index].sample(len(CPGs_in_index), random_state = 14)
    
#     return CPGs_subset, non_CPGs_subset




# Load CPGs and CPG-like metadata
CPG_metadata = pd.read_csv(snakemake.input[1], sep = '\t', header = 0)


# Load frequencies table
df = pd.read_csv(input_file, sep = '\t', header = 0, index_col = 0)


# Get CPG subset
if 'CPGs' in output_file:
    CPGs_freqs_subset = CPGs_subset(df, CPG_metadata)
    CPGs_freqs_subset.to_csv( path_or_buf = output_file, sep = '\t', index =  True)

# Get nonCPG subset
if 'nonCPGs' in output_file:
    nonCPGs_freqs_subset = nonCPGs_subset(df, CPG_metadata)
    nonCPGs_freqs_subset.to_csv( path_or_buf = output_file, sep = '\t', index =  True)