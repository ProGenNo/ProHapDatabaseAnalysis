configfile: "config.yaml"

ENZYMES = ["Trypsin", "Lys-C", "Lys-N", "Glu-C", "Asp-N"]
POPULATIONS = ["EUR", "EAS", 'SAS', 'AMR', 'AFR']

rule all:
    input:
        pept=config['final_peptide_list'],
        discoverable_vars=config['discoverable_variant_list'],
        haplo_vars=expand('{proxy}', proxy=[config['possible_variant_list']] if len(config["haplo_db_table"]) > 0 else []),
        populations_proteome_coverage=expand("results/pep_coverage_{popul}.py", popul=POPULATIONS)

rule list_all_possible_variants:
    input:
        config['haplo_db_table']
    output:
        v=config['possible_variant_list']
    conda: "envs/main_env.yaml"
    shell:
        "python3 src/haplo_extract_all_vars.py -hap_tsv {input} -o {output.v} "

rule digest_proteins:
    input:
        config['full_fasta_file']
    output:
        temp('results/peptide_list_{enz}.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "mkdir -p results; "
        "python3 src/create_peptide_list.py -i {input} -enz {wildcards.enz} -o {output}"

rule merge_peptide_lists:
    input:
        expand('results/peptide_list_{enz}.tsv', enz=ENZYMES)
    output:
        "results/peptide_list_full.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand('results/peptide_list_{enz}.tsv', enz=ENZYMES))
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule peptides_annotate_variation:
    input:
        peptides="results/peptide_list_full.tsv",
        var_db=expand('{proxy}', proxy=[config['var_db_table']] if len(config["var_db_table"]) > 0 else []),
        haplo_db=expand('{proxy}', proxy=[config['haplo_db_table']] if len(config["haplo_db_table"]) > 0 else []),
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['full_fasta_file'],
        ref_fasta=config['reference_fasta']
    output:
        config['final_peptide_list']
    params:
        variant_prefix=config['variant_protein_prefix'],
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/variation_annot_errors.log"
    threads: config['max_cores']
    #conda: "envs/main_env.yaml"
    run:
        if ((len(config["var_db_table"]) > 0) and (len(config["haplo_db_table"]) > 0)):
            shell("python3 src/peptides_annotate_variation.py -i {input.peptides} \
                        -var_tsv {input.var_db} -hap_tsv {input.haplo_db} \
                        -var_prefix {params.variant_prefix} -hap_prefix {params.haplotype_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")
    
        elif (len(config["var_db_table"]) > 0):
            shell("python3 src/peptides_annotate_variation.py -i {input.peptides} \
                        -var_tsv {input.var_db} -var_prefix {params.variant_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")

        elif (len(config["haplo_db_table"]) > 0):
            shell("python3 src/peptides_annotate_variation.py -i {input.peptides} \
                        -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")

rule get_discoverable_variants:
	input:
		pep=config['final_peptide_list'],
		hap=config['haplo_db_table']
	output:
		config['discoverable_variant_list']
	conda: "envs/main_env.yaml"
	shell:
		"python src/get_discoverable_variants.py -i {input.pep} -hap_tsv {input.hap} -o {output}; "

rule get_coverage:
    input:
        pep=config['final_peptide_list'],
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file=config['full_fasta_file'],
        ref_fasta=config['reference_fasta']
    output:
        "results/pep_coverage.py"
    params:
        max_cores=config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_stats_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -tr_id {input_tr_ids} -o {output}"

rule filter_to_population:
    input:
        hap=config['haplo_db_table'],
        fasta=config['full_fasta_file']
    output:
        hap="data/haplotypes_{popul}.tsv",
        fasta="data/proteindb_{popul}.tsv"
    params:
        max_cores=config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/filter_haplotypes.py -f {input.fasta} -hap_tsv {input.hap} -pop {wildcards.pop} -t {params.max_cores} -output_tsv {output.hap} -output_fasta {output.fasta}"

rule digest_proteins_pop:
    input:
        "data/proteindb_{popul}.tsv"
    output:
        temp('results/peptide_list_{enz}_{popul}.tsv')
    conda: "envs/main_env.yaml"
    shell:
        "mkdir -p results; "
        "python3 src/create_peptide_list.py -i {input} -enz {wildcards.enz} -o {output}"
    
rule merge_peptide_lists_pop:
    input:
        expand('results/peptide_list_{enz}_{{popul}}.tsv', enz=ENZYMES)
    output:
        "results/peptide_list_full_{popul}.tsv"
    conda: "envs/main_env.yaml"
    params:
        input_file_list = ','.join(expand('results/peptide_list_{enz}_{{popul}}.tsv', enz=ENZYMES))
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule peptides_annotate_variation_pop:
    input:
        peptides="results/peptide_list_full_{popul}.tsv",
        haplo_db="data/haplotypes_{popul}.tsv",
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file="data/proteindb_{popul}.tsv",
        ref_fasta=config['reference_fasta']
    output:
        "results/peptide_list_{popul}.csv"
    params:
        variant_prefix=config['variant_protein_prefix'],
        haplotype_prefix=config['haplotype_protein_prefix'],
        max_cores=config['max_cores'],
        log_file="results/variation_annot_errors.log"
    threads: config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/peptides_annotate_variation.py -i {input.peptides} -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; "

rule get_coverage_pop:
    input:
        pep="results/peptide_list_{popul}.csv",
        tr_ids='data/protein_transcript_ids_110.csv',
        gene_ids='data/gene_transcript_ids_110.csv',
        fasta_file="data/proteindb_{popul}.tsv",
        ref_fasta=config['reference_fasta']
    output:
        "results/pep_coverage_{popul}.py"
    params:
        max_cores=config['max_cores']
    conda: "envs/main_env.yaml"
    shell:
        "python src/get_peptide_stats_parallel.py -i {input.pep} -f {input.fasta_file} -t {params.max_cores} -ref_fa {input.ref_fasta} -g_id {input.gene_ids} -tr_id {input_tr_ids} -o {output}"

    