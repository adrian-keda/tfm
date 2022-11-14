configfile: "config.yaml"



#Define parameters
chr = config["chromosomes"]
thr = config["threshold"]
patho = config["pathogenic"]
subset = config["gene_subset"]
out = config["outdir"]


rule all:
    input:
        expand(f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.{{subset}}_frequencies.tsv", patho = patho, subset = subset)



rule s1_clean_information:
    input:
        f"Chromosomes/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv",
        config["files"]["ancestry"],
        config["files"]["pop_count"]
    output:
        f"{out}clean_info/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.{{chr}}.tsv"
    script:
        "scripts/0_clean_popcancer_freqs.py"



rule s1_5_coding_variants_filter:
    input:
        f"{out}clean_info/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.{{chr}}.tsv"
    output:
        f"{out}coding/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.{{chr}}.tsv"
    script:
        "scripts/0.5_filter_coding_variants.py"



rule s2_rare_variants_filter:
    input:
        f"{out}coding/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.{{chr}}.tsv"
    output:
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.{{chr}}.tsv"
    params: thr
    script:
        "scripts/1_filter_rare_variants.py"



rule s3_get_header:
    input:
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.chr21.tsv"
    output:
        f"{out}header.tsv"
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
        f"{out}rare_vars/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.clean.popcancer.freqs.canonic.conseq.coding.rare_{thr}.all_chr.header.tsv"
    output:
        f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.tsv"
    wildcard_constraints:
        patho = "([A-Z]|[a-z]){6,7}"
    script:
        "scripts/2_filter_pathogenic_variants.py"



rule s7_table_gene_vs_cancer:
    input:
        f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.tsv",
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
        f"{out}pathogenic_vars_fr_eq_{thr}/rare_vars_all_chr.pathogenic_{{patho}}.{{subset}}_frequencies.tsv"
    script:
        "scripts/4_CPGs_nonCPGs_frequencies.py"
