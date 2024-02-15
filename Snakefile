configfile: "config.yaml"

ENZYMES = ["Trypsin", "Lys-C", "Lys-N", "Glu-C", "Asp-N", "Chymotrypsin"]
POPULATIONS = ["EUR", "EAS", 'SAS', 'AMR', 'AFR']

rule all:
    input:
        pept="results/peptide_list_ALL.tsv",
        discoverable_vars="results/all_discoverable_variants.tsv",
        var_stats="results/variant_stats.txt",
        haplo_vars="results/all_included_variants.csv",
        proteome_coverage="results/peptide_coverage_stats.tsv",
        freq_coverage="results/peptide_coverage_freq_log.tsv",
        freq_per_transcript="results/transcript_freqs_by_superpop.tsv"

rule list_all_possible_variants:
    input:
        config['haplo_db_table']
    output:
        v="results/all_included_variants.csv"
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/haplo_extract_all_vars.py -hap_tsv {input} -o {output.v} "

rule digest_proteins:
    input:
        config['full_fasta_file']
    output:
        temp('results/peptide_list_{enz}.tsv')
    params:
        missed_cl=lambda wildcards: 4 if (wildcards.enz == 'Chymotrypsin') else 2
    conda: "envs/main_env.yaml"
    shell:
        "mkdir -p results; "
        "python3 src/create_peptide_list.py -i {input} -enz {wildcards.enz} -o {output} -m {params.missed_cl}"

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
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_list_ALL.tsv"
    params:
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/variation_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule get_discoverable_variants:
	input:
		pep="results/peptide_list_ALL.tsv",
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
        pep=config['final_peptide_list'],
        tr_ids='data/protein_transcript_ids_110.csv',
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
        "python src/get_peptide_stats_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -tr_id {input.tr_ids} -o {output}"

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
        fasta_file="data/proteindb_{popul}.fa",
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
        "python src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule add_peptide_frequency:
    input:
        pep="results/peptide_temp_list_{popul}.csv",
        haplo_db="data/haplotypes_{popul}.tsv"
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
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv'
    output:
        "results/pep_freq_coverage_{popul}.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_coverage_freq.py -i {input.pep} -t {params.max_cores} -g_id {input.gene_ids} -tr_id {input.tr_ids} -o {output}"

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
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file="data/proteindb_{popul}.fa",
        ref_fasta=config['reference_fasta']
    output:
        "results/pep_coverage_{popul}.tsv"
    params:
        max_cores=config['max_cores']
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_stats_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -tr_id {input.tr_ids} -o {output}"

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
        input_file_list = ','.join(expand(config['popul_haplo_table'], popul=POPULATIONS)),
        max_cores=config['max_cores']
    threads: config['max_cores']
    shell:
        "python src/extract_population_frequencies.py -i {params.input_file_list} -transcripts {input.tr_id} -s {input.samples} -t {params.max_cores} -o_superpop {output}"
    
    