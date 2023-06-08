configfile: "config.yaml"

ENZYMES = ["Trypsin", "Lys-C", "Lys-N", "Glu-C", "Asp-N", "Chymotrypsin"]

rule all:
    input:
        pept=config['final_peptide_list']

rule digest_proteins:
    input:
        config['full_fasta_file']
    output:
        temp('results/peptide_list_{enz}.tsv')
    shell:
        "mkdir -p results; "
        "python3 src/analysis/create_peptide_list.py -i {input} -enz {wildcards.enz} -o {output}"

rule merge_peptide_lists:
    input:
        expand('results/peptide_list_{enz}.tsv', enz=ENZYMES)
    output:
        temp("results/peptide_list_full.tsv")
    params:
        input_file_list = ','.join(expand('results/peptide_list_{enz}.tsv', enz=ENZYMES))
    shell:
        "python3 src/merge_tables.py -i {params.input_file_list} -o {output}"

rule annotate_variation:
    input:
        peptides="results/peptide_list_full.tsv",
        var_db=expand('{proxy}', proxy=[config['var_db_table']] if len(config["var_db_table"]) > 0 else []),
        haplo_db=expand('{proxy}', proxy=[config['haplo_db_table']] if len(config["haplo_db_table"]) > 0 else []),
        tr_ids='data/protein_transcript_ids_108.csv',
        gene_ids='data/gene_transcript_ids_108.csv',
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
    run:
        if ((len(config["var_db_table"]) > 0) and (len(config["haplo_db_table"]) > 0)):
            shell("python3 src/results_postprocess_newDB/psm_annotate_variation_newDB.py -i {input.peptides} \
                        -var_tsv {input.var_db} -hap_tsv {input.haplo_db} \
                        -var_prefix {params.variant_prefix} -hap_prefix {params.haplotype_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")
    
        elif (len(config["var_db_table"]) > 0):
            shell("python3 src/results_postprocess_newDB/psm_annotate_variation_newDB.py -i {input.peptides} \
                        -var_tsv {input.var_db} -var_prefix {params.variant_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")

        elif (len(config["haplo_db_table"]) > 0):
            shell("python3 src/results_postprocess_newDB/psm_annotate_variation_newDB.py -i {input.peptides} \
                        -hap_tsv {input.haplo_db} -hap_prefix {params.haplotype_prefix} \
                        -log {params.log_file} -tr_id {input.tr_ids} -g_id {input.gene_ids} -f {input.fasta_file} -ref_fa {input.ref_fasta} -t {params.max_cores} -o {output}; ")