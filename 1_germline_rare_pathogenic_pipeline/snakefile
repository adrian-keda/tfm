configfile: "config.yaml"



#Define parameters
chr = config["chromosomes"]
thr = config["threshold"]
patho = config["pathogenic"]
subset = config["gene_subset"]
out = config["outdir"]


rule all:
    input:
        expand(f"{out}important_things_{thr}/pathogenic_vars/rare_vars_all_chr.pathogenic_{{patho}}_{thr}.tsv", patho = patho),
        expand(f"{out}important_things_{thr}/tables/rare_vars_all_chr.pathogenic_{{patho}}.{{subset}}_freqs_{thr}.tsv", patho = patho, subset = subset),
        f"{out}important_things_{thr}/log_plots/general_filtering_log_{thr}.tsv",
        expand(f"{out}important_things_{thr}/log_plots/log_plot_{{patho}}_{thr}.svg", patho = patho)



rule s1_clean_information:
    input:
        input = f"Chromosomes/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv",
        ancestry = config["files"]["ancestry"],
        pop_count = config["files"]["pop_count"]
    output:
        output = temp(f"{out}clean_info/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.{{chr}}.tsv"),
        vars_filtered_log = temp(f"{out}variants_filtered_log/clean/{{chr}}.tsv")
    script:
        "scripts/0_clean_popcancer_freqs.py"



rule s1_5_coding_variants_filter:
    input:
        input = f"{out}clean_info/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.{{chr}}.tsv",
        vars_filtered_log = f"{out}variants_filtered_log/clean/{{chr}}.tsv"
    output:
        output = f"{out}coding/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.{{chr}}.tsv",
        vars_filtered_log = f"{out}variants_filtered_log/coding/{{chr}}.tsv"
    script:
        "scripts/0.5_filter_coding_variants.py"


rule s2_rare_variants_filter:
    input:
        input = f"{out}coding/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.{{chr}}.tsv",
        vars_filtered_log = f"{out}variants_filtered_log/coding/{{chr}}.tsv"
    output:
        output = f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.{{chr}}.tsv",
        vars_filtered_log = f"{out}variants_filtered_log/rare_{thr}/{{chr}}.tsv"
    params: thr
    script:
        "scripts/1_filter_rare_variants.py"



rule s3_get_header:
    input:
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.chr21.tsv"
    output:
        temp(f"{out}header.tsv")
    shell:
        "head -n 1 {input} > {output}"



rule s3_5_remove_header:
    input:
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.{{chr}}.tsv"
    output:
        temp(f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.{{chr}}.noheader.tsv")
    shell:
        "tail -n +2 {input} > {output}"



rule s4_concat_files:
    input:
        expand(f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.{{chr}}.noheader.tsv", chr = chr)
    output:
        temp(f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.all_chr.tsv")
    shell:
        "cat {input} >> {output}"



rule s5_add_header:
    input:
        all = f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.all_chr.tsv",
        header = f"{out}header.tsv"
    output:
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.all_chr.header.tsv"
    shell:
        "cat {input.header} {input.all} >> {output}"



rule s6_pathogenic_variants_filter:
    input:
        input = f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.all_chr.header.tsv",
    output:
        output = f"{out}important_things_{thr}/pathogenic_vars/rare_vars_all_chr.pathogenic_{{patho}}_{thr}.tsv",
        vars_filtered_log = f"{out}variants_filtered_log/{{patho}}_{thr}.tsv"
    wildcard_constraints:
        patho = "([A-Z]|[a-z]){6,7}"
    script:
        "scripts/2_filter_pathogenic_variants.py"



rule s7_table_gene_vs_cancer:
    input:
        f"{out}important_things_{thr}/pathogenic_vars/rare_vars_all_chr.pathogenic_{{patho}}_{thr}.tsv",
        config["files"]["cancer_count"]
    output:
        f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.frequencies.tsv"
    script:
        "scripts/3_table_gene_vs_cancer.py"



rule s8_substract_frequencies:
    input:
        f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.frequencies.tsv",
        config["files"]["CPG_CPGlike_list"]
    output:
        f"{out}important_things_{thr}/tables/rare_vars_all_chr.pathogenic_{{patho}}.{{subset}}_freqs_{thr}.tsv"
    script:
        "scripts/4_CPGs_nonCPGs_frequencies.py"



rule s9_merging_logs:
    input:
        chr_logs = expand(f"{out}variants_filtered_log/rare_{thr}/{{chr}}.tsv", chr = chr),
        patho_logs = expand(f"{out}variants_filtered_log/{{patho}}_{thr}.tsv", patho = patho)
    output:
        output_log = f"{out}important_things_{thr}/log_plots/general_filtering_log_{thr}.tsv"
    script:
        "scripts/5_merging_logs.py"


rule s10_draw_plots:
    input:
        pipeline_log = f"{out}important_things_{thr}/log_plots/general_filtering_log_{thr}.tsv"
    output:
        expand(f"{out}important_things_{thr}/log_plots/log_plot_{{patho}}_{thr}.svg", patho = patho)
    script:
        "scripts/6_variants_plot.py"