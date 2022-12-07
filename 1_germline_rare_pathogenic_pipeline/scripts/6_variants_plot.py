import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




# Reading log file
df = pd.read_csv(snakemake.input["pipeline_log"], sep = '\t', header = 0)


# # LETS DO SOME FANCY PLOTS
# def number_plot(df):

#     plt.figure(figsize = (9,7))
#     # Define values to plot
#     x = range(1,6)
#     y1 = df['Input'].tolist()[0:5]
#     y2 = df['Output'].tolist()[0:5]

#     # Plot itself
#     plt.plot(x, y1, marker = 'o', color = 'blue')
#     plt.plot(x, y2, marker = '^', color = 'red')

#     # Text
#     plt.xlabel('Filtering step', fontweight = 'bold', fontsize = 'large')
#     plt.ylabel('Number of variants', fontweight = 'bold', fontsize = 'large', rotation = 'vertical')
#     plt.xticks(ticks = x, labels = df['Filter_ID'].tolist()[0:5], fontsize = 'small', rotation = 15)

#     plt.legend(labels = ['Input variants', 'Output variants'])
#     plt.title('Number of variants per filtering step', fontsize = 17, fontweight = 'bold')

#     plt.savefig(snakemake.output["number_plot"])




# def percentage_plot(df):
    
#     plt.figure(figsize = (9,7))
#     # Define values to plot
#     x = range(1,8)
#     input_vars = df['Input'].tolist()[0]
#     percentage = df['Difference'] / input_vars
    
#     labels = df['Filter_ID'].tolist()

#     # Plot itself
#     plt.ylim(0, max(percentage))
#     plt.bar(x, percentage, width = 0.55, align = 'center', color = '#8c8c8c', edgecolor = 'black')

#     # Text options
#     plt.xlabel('Filtering step', fontweight = 'bold', fontsize = 'large')
#     plt.ylabel('Percentage', fontweight = 'bold', fontsize = 'large', rotation = 'vertical')
#     for i in x:
#         plt.text(x = i, y = 0.05, s = round(percentage [i-1], 4), fontweight = 'bold', horizontalalignment = 'center')
#     plt.xticks(ticks = x, labels = labels, fontsize = 'small', rotation = 15)
#     plt.yticks(ticks = [0.1, 0.2, 0.3, 0.4, 0.5, round(max(percentage), 2)])
#     plt.title('Percentage of input variants filtered per step', fontsize = 17, fontweight = 'bold')

#     plt.savefig(snakemake.output["percentage_plot"])




def MASTER_PLOTTING_FUNCTION_OMEGALUL_XDDD (df, output_name):

    if 'DelMis' in output_name:
        patho = 'DelMis'
    else:
        patho = 'ClinVar'

    plt.figure(figsize = (10,7), dpi = 100)

    # Define values to plot
    input_vars = df['Input'].tolist()[0]
    output_vars = df['Output'].tolist()[0:5]
    patho_vars = df[df['Filter_ID'] == patho]['Output'].tolist()
    data = [input_vars] + output_vars + patho_vars

    data_formated = ["{:.2e}".format(i) for i in data]
    percentages = [round((i / input_vars) * 100, 2) for i in data]

    labels = ['Input vars'] + df['Filter_ID'].tolist()[0:5] + [patho]
    x = np.arange(1, len(labels) + 1).tolist()

    # Plot
    plt.plot(x, percentages, marker = 'o', color = '#d11141', markersize = 7, linewidth = 1.5, zorder = 4, clip_on = False)
    plt.ylim (0, 100)

    # Text options
    #pos = 0
    for i, txt in enumerate (data_formated):
        plt.annotate(f'{str(txt)}\n{percentages[i]}%', (x[i], percentages[i]), fontweight = 'bold', clip_on = False, zorder = 5)
    plt.xlabel('Filtering step', fontweight = 'bold', fontsize = 15)
    plt.ylabel('% of variants left', fontweight = 'bold', fontsize = 15, rotation = 'vertical')
    plt.xticks(ticks = x, labels = labels, fontsize = 9, rotation = 15)
    plt.title('Variants left per filtering step', y = 1.05, fontsize = 25, fontweight = 'bold')

    plt.savefig(output_name, format = 'svg')


# number_plot(df)
# percentage_plot(df)
for name in snakemake.output:
    MASTER_PLOTTING_FUNCTION_OMEGALUL_XDDD(df, name)