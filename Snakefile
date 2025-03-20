configfile: "config.yaml"

ENZYMES = ["Trypsin", "Lys-C", "Lys-N", "Glu-C", "Asp-N", "Chymotrypsin"]
POPULATIONS = ["EUR", "EAS", 'SAS', 'AMR', 'AFR']

rule all:
    input:
        pept="results/peptide_list_ALL.csv",
        discoverable_vars="results/all_discoverable_variants.tsv",
        var_stats="results/variant_stats.txt",
        haplo_vars="results/all_included_variants.csv",
        proteome_coverage="results/peptide_coverage_stats.tsv",
        freq_coverage="results/peptide_coverage_freq_log.tsv",
        indiv_coverage="results/peptide_indiv_coverage.tsv",
        freq_per_transcript="results/transcript_freqs_by_superpop.tsv",
        uniprot_comparison="results/uniprot_comparison_stats.tsv"

rule list_all_possible_variants:
    input:
        config['haplo_db_table']
    output:
        "results/all_included_variants.csv"
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/haplo_extract_all_vars.py -hap_tsv {input} -o {output} "

rule digest_proteins:
    input:
        fasta=config['full_fasta_file']
    output:
        temp('results/peptide_list_{enz}.tsv')
    params:
        missed_cl=lambda wildcards: 4 if (wildcards.enz == 'Chymotrypsin') else 2
    conda: "envs/main_env.yaml"
    shell:
        "mkdir -p results; "
        "python3 src/create_peptide_list.py -i {input.fasta} -enz {wildcards.enz} -o {output} -m {params.missed_cl}"

rule merge_peptide_lists:
    input:
        expand('results/peptide_list_{enz}.tsv', enz=ENZYMES)
    output:
        temp("results/peptide_list_full.tsv")
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand('results/peptide_list_{enz}.tsv', enz=ENZYMES))
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule peptides_annotate_variation:
    input:
        peptides="results/peptide_list_full.tsv",
        haplo_db=config["haplo_db_table"],
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['full_fasta_file'],
        fasta_header=config['full_fasta_header'],
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_list_ALL.csv"
    params:
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/variation_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -fh {input.fasta_header} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule get_discoverable_variants:
	input:
		pep="results/peptide_list_ALL.csv",
		hap=config['haplo_db_table']
	output:
		"results/all_discoverable_variants.tsv"
	conda: "envs/main_env.yaml"
	shell:
		"python src/get_discoverable_variants.py -i {input.pep} -hap_tsv {input.hap} -o {output}; "

rule collect_variant_stats:
    input:
        all="results/all_included_variants.csv",
        discoverable="results/all_discoverable_variants.tsv"
    output:
        "results/variant_stats.txt"
    shell:
        "python src/get_variant_counts.py -i_all {input.all} -i_disc {input.discoverable} > {output}"

rule get_coverage:
    input:
        pep="results/peptide_list_ALL.csv",
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['full_fasta_file'],
        ref_fasta=config['reference_fasta']
    output:
        "results/pep_coverage.tsv"
    params:
        max_cores=config['max_cores']
    conda: "envs/main_env.yaml"
    threads: config['max_cores']
    shell:
        "python src/get_peptide_coverage_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -o {output}"

rule digest_proteins_pop:
    input:
        config['population_fasta_file']
    output:
        temp('results/pop_peptide_list_{enz}_{popul}.tsv')
    params:
        missed_cl=lambda wildcards: 4 if (wildcards.enz == 'Chymotrypsin') else 2
    conda: "envs/main_env.yaml"
    shell:
        "mkdir -p results; "
        "python3 src/create_peptide_list.py -i {input} -enz {wildcards.enz} -o {output} -m {params.missed_cl}"
    
rule merge_peptide_lists_pop:
    input:
        expand('results/pop_peptide_list_{enz}_{{popul}}.tsv', enz=ENZYMES)
    output:
        "results/pop_{popul}_peptide_list_full.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand('results/pop_peptide_list_{enz}_{{popul}}.tsv', enz=ENZYMES))
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule peptides_annotate_variation_pop:
    input:
        peptides="results/pop_{popul}_peptide_list_full.tsv",
        haplo_db=config['popul_haplo_table'],
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['population_fasta_file'],
        fasta_header=config['population_fasta_header'],
        ref_fasta=config['reference_fasta']
    output:
        temp("results/peptide_temp_list_{popul}.csv")
    params:
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/variation_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -fh {input.fasta_header} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule add_peptide_frequency:
    input:
        pep="results/peptide_temp_list_{popul}.csv",
        haplo_db=config['popul_haplo_table']
    output:
        "results/peptide_list_{popul}.csv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/peptides_add_frequency.py -i {input.pep} -hap_tsv {input.haplo_db} -pop {wildcards.popul} -t {params.max_cores} -o {output}"

rule get_coverage_freq_pop:
    input:
        pep="results/peptide_list_{popul}.csv",
        gene_ids='data/gene_transcript_ids_110.csv'
    output:
        "results/pep_freq_coverage_{popul}.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_coverage_freq.py -i {input.pep} -t {params.max_cores} -g_id {input.gene_ids} -o {output}"

rule collect_coverage_freq_stats:
    input:
        pop_coverage=expand("results/pep_freq_coverage_{popul}.tsv", popul=POPULATIONS),
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_coverage_freq_log.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand("results/pep_freq_coverage_{popul}.tsv", popul=POPULATIONS)),
        populations = ','.join(POPULATIONS)
    shell:
        "python src/get_peptide_stats_aggregated_freq.py -i {params.input_file_list} -pop {params.populations} -ref_fa {input.ref_fasta} -o {output}"

rule get_coverage_pop:
    input:
        pep="results/peptide_list_{popul}.csv",
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['population_fasta_file'],
        ref_fasta=config['reference_fasta']
    output:
        "results/pep_coverage_{popul}.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_coverage_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -o {output}"

rule collect_coverage_stats:
    input:
        pop_coverage=expand("results/pep_coverage_{popul}.tsv", popul=POPULATIONS),
        all_coverage="results/pep_coverage.tsv",
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_coverage_stats.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(["results/pep_coverage.tsv"] + expand("results/pep_coverage_{popul}.tsv", popul=POPULATIONS)),
        populations = ','.join(['all'] + POPULATIONS)
    shell:
        "python src/get_peptide_stats_aggregated.py -i {params.input_file_list} -pop {params.populations} -ref_fa {input.ref_fasta} -o {output}"
    
rule get_frequencies_per_transcript:
    input:
        all_haplotypes=expand(config['popul_haplo_table'], popul=POPULATIONS),
        tr_id='data/protein_transcript_ids_110.csv',
        samples='igsr_samples_filtered.tsv'
    output:
        "results/transcript_freqs_by_superpop.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_filename = lambda wildcards: config['popul_haplo_table'],
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python src/extract_population_frequencies.py -i {params.input_filename} -transcripts {input.tr_id} -s {input.samples} -t {params.max_cores} -o_superpop {output}"
    
rule get_individual_coverage:
    input:    
        pep="results/peptide_list_{popul}.csv",
        haplo_db=config['haplo_samples_table'],
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['population_fasta_file'],
        ref_fasta=config['reference_fasta'],
        samples='igsr_samples_filtered.tsv'
    output:
        "results/pep_indiv_coverage_{popul}.tsv"
    conda: "envs/main_env.yaml"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python src/get_peptide_coverage_indiv.py -i {input.pep} -t {params.max_cores} -g_id {input.gene_ids} -o {output} -s {input.samples} -hap_tsv {input.haplo_db} -pop {wildcards.popul}"

rule collect_coverage_stats_indiv:
    input:
        pop_coverage=expand("results/pep_indiv_coverage_{popul}.tsv", popul=POPULATIONS),
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_indiv_coverage.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand("results/pep_indiv_coverage_{popul}.tsv", popul=POPULATIONS)),
    shell:
        "python src/get_peptide_stats_aggregated_indiv.py -i {params.input_file_list} -ref_fa {input.ref_fasta} -o {output} "
    
rule unzip_swissprot:
    input:
        'data/uniprot/uniprotkb_proteome_UP000005640_AND_revi_2024_06_25.fasta.gz'
    output:
        temp('data/uniprot/uniprotkb_proteome_UP000005640_AND_revi_2024_06_25.fasta')
    shell:
        'gunzip {input}'

rule unzip_uniprot:
    input:
        'data/uniprot/uniprotkb_proteome_UP000005640_2024_06_25_incl_unreviewed.fasta.gz'
    output:
        temp('data/uniprot/uniprotkb_proteome_UP000005640_2024_06_25_incl_unreviewed.fasta')
    shell:
        'gunzip {input}'

rule digest_proteins_swissprot:
    input:
        'data/uniprot/uniprotkb_proteome_UP000005640_AND_revi_2024_06_25.fasta'
    output:
        temp('data/uniprot/swissprot_trypsin.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/create_peptide_list.py -i {input} -enz Trypsin -o {output} -m 2"

rule digest_proteins_uniprot:
    input:
        'data/uniprot/uniprotkb_proteome_UP000005640_2024_06_25_incl_unreviewed.fasta'
    output:
        temp('data/uniprot/uniprot_trypsin.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/create_peptide_list.py -i {input} -enz Trypsin -o {output} -m 2"

rule digest_proteins_hrc:
    input:
        config['hrc_fasta_file']
    output:
        temp('data/HRC/hrc_peptides_trypsin.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/create_peptide_list.py -i {input} -enz Trypsin -o {output} -m 2"

rule annotate_peptides_hrc:
    input:
        peptides="data/HRC/hrc_peptides_trypsin.tsv",
        haplo_db=config["hrc_haplo_table"],
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['hrc_fasta_file'],
        fasta_header=config['hrc_fasta_header'],
        ref_fasta=config['reference_fasta']
    output:
        temp('data/HRC/hrc_peptides_trypsin_annot.tsv')
    params:
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/hrc_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -fh {input.fasta_header} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule digest_proteins_pangenome:
    input:
        config['pangen_fasta_file']
    output:
        temp('data/pangenome_peptides_trypsin.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/create_peptide_list.py -i {input} -enz Trypsin -o {output} -m 2"
    
rule annotate_peptides_pangenome:
    input:
        peptides='data/pangenome_peptides_trypsin.tsv',
        haplo_db=config["pangen_haplo_table"],
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['pangen_fasta_file'],
        fasta_header=config['pangen_fasta_header'],
        ref_fasta=config['reference_fasta']
    output:
        temp('data/pangenome_peptides_trypsin_annot.tsv')
    params:
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/hrc_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -fh {input.fasta_header} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule compare_uniprot_hrc:
    input:
        sp='data/uniprot/swissprot_trypsin.tsv',
        up='data/uniprot/uniprot_trypsin.tsv',
        pep='data/HRC/hrc_peptides_trypsin_annot.tsv',
    output:
        temp('results/uniprot_comparison_peptides_hrc.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python src/compare_with_uniprot.py -i {input.pep} -sp {input.sp} -iso {input.up} -o {output}"

rule compare_uniprot_pangenome:
    input:
        sp='data/uniprot/swissprot_trypsin.tsv',
        up='data/uniprot/uniprot_trypsin.tsv',
        pep='data/pangenome_peptides_trypsin_annot.tsv',
    output:
        temp('results/uniprot_comparison_peptides_pangenome.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python src/compare_with_uniprot.py -i {input.pep} -sp {input.sp} -iso {input.up} -o {output}"

rule compare_uniprot_pop:
    input:
        sp='data/uniprot/swissprot_trypsin.tsv',
        up='data/uniprot/uniprot_trypsin.tsv',
        pop_pep='results/peptide_list_{popul}.csv',
    output:
        temp('results/uniprot_comparison_peptides_{popul}.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python src/compare_with_uniprot.py -i {input.pop_pep} -sp {input.sp} -iso {input.up} -o {output}"
    
rule compare_uniprot_1kgp:
    input:
        sp='data/uniprot/swissprot_trypsin.tsv',
        up='data/uniprot/uniprot_trypsin.tsv',
        pep='results/peptide_list_ALL.csv',
    output:
        temp('results/uniprot_comparison_peptides_1kgp.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "python src/compare_with_uniprot.py -i {input.pep} -sp {input.sp} -iso {input.up} -o {output}"

rule collect_uniprot_comparison:
    input:
        pop_pep=expand("results/uniprot_comparison_peptides_{popul}.tsv", popul=POPULATIONS),
        all_pep='results/uniprot_comparison_peptides_1kgp.tsv',
        hrc_pep='results/uniprot_comparison_peptides_hrc.tsv',
        pangen_pep="results/uniprot_comparison_peptides_pangenome.tsv"
    output:
        "results/uniprot_comparison_stats.tsv"    
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(["results/uniprot_comparison_peptides_1kgp.tsv", "results/uniprot_comparison_peptides_pangenome.tsv", 'results/uniprot_comparison_peptides_hrc.tsv'] + expand("results/uniprot_comparison_peptides_{popul}.tsv", popul=POPULATIONS)),
        populations = ','.join(['ALL', 'Pangenome_ALL', 'HRC'] + POPULATIONS)
    shell:
        "python src/get_uniprot_comparison_aggregated.py -i {params.input_file_list} -pop {params.populations} -o {output}"
    
rule check_IL_ambiguities:
    input:
        sp='data/uniprot/swissprot_trypsin.tsv',
        up='data/uniprot/uniprot_trypsin.tsv',
        pep='results/peptide_list_ALL.csv',
    output:
        "results/IL_ambiguities_ProHap_SP_UP_Trypsin.csv"
    conda: "envs/main_env.yaml"
    shell:
        "python src/check_IL_ambiguities.py -i {input.pep} -sp {input.sp} -iso {input.pt} -o {output}"
