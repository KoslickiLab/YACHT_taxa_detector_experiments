configfile: "config.yaml"

wildcard_constraints:
    dataset="\d"

rule download_sample:
    output:
        "sample_files/sample_file{num}.tar.gz"
    shell:
        "mkdir -p sample_files; wget -O ./{input} https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_{wildcards.num}_reads.tar.gz"

rule unzip_sample_step1:
    input:
        "sample_files/sample_file{num}.tar.gz"
    output:
        "sample_folders/sample{num}/simulation_short_read/2018.08.15_09.49.32_sample_{num}/reads/anonymous_reads.fq.gz",
        "sample_folders/sample{num}/simulation_short_read/2018.08.15_09.49.32_sample_{num}/reads/reads_mapping.tsv.gz"
    shell:
        "mkdir -p sample_folders/sample{wildcards.num}; tar -xvzf ./{input} -C ./sample_folders/sample{wildcards.num}"

rule unzip_sample_step2:
    input:
        "sample_folders/sample{num}/simulation_short_read/2018.08.15_09.49.32_sample_{num}/reads/anonymous_reads.fq.gz",
        "sample_folders/sample{num}/simulation_short_read/2018.08.15_09.49.32_sample_{num}/reads/reads_mapping.tsv.gz"
    output:
        "sample_folders/sample{num}/anonymous_reads.fq",
        "sample_folders/sample{num}/reads_mapping.tsv"
    shell:
        "cp {input[0]} {output[0]}.gz; cp {input[1]} {output[1]}.gz; gzip -dv {output[0]}.gz; gzip -dv {output[1]}.gz"

scaled = ("100" if config["use_ncbi_database"] else "1000")

rule sketch_sample:
    input:
        "sample_folders/sample{num}/anonymous_reads.fq"
    output:
        "sample{num}.sig.zip"
    conda:
        "yacht_env"
    shell:
        "yacht sketch sample --infile {input} --kmer 31 --scaled {scaled} --outfile {output}"

rule download_gtdb_ref:
    output:
        "ref_gtdb/gtdb-rs214-reps.k31_0.9995_pretrained/gtdb-rs214-reps.k31_0.9995_config.json"
    conda:
        "yacht_env"
    shell:
        "yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 31 --ani_thresh 0.9995 --outfolder ref_gtdb"

ref_db = ("ref_ncbi/customized_ncbi_ani_thresh_0.95_config.json" if config["use_ncbi_database"] else "ref_gtdb/gtdb-rs214-reps.k31_0.9995_pretrained/gtdb-rs214-reps.k31_0.9995_config.json")

rule run_yacht:
    input:
        "sample{num}.sig.zip",
        ref_db
    output:
        "result_sample{num}.xlsx"
    conda:
        "yacht_env"
    shell:
        "yacht run --json '{input[1]}' --sample_file '{input[0]}' --num_threads 64 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out {output}"

genome_to_taxid_script = ("create_genome_to_taxid_ncbi.py" if config["use_ncbi_database"] else "create_genome_to_taxid_gtdb.py")

rule create_genome_to_taxid:
    input:
        "result_sample{num}.xlsx"
    output:
        "genome_to_taxid_sample{num}.tsv"
    conda:
        "cami_env.yaml"
    script:
        genome_to_taxid_script

rule convert_yacht_output:
    input:
        "result_sample{num}.xlsx",
        "genome_to_taxid_sample{num}.tsv"
    output:
        "cami_result_sample{num}.cami"
    conda:
        "yacht_env"
    shell:
        "yacht convert --yacht_output '{input[0]}' --sheet_name 'min_coverage0.1' --genome_to_taxid '{input[1]}' --mode 'cami' --sample_name 'marmgCAMI2_short_read_sample_0' --outfile_prefix 'cami_result' --outdir ./"

rule download_gt:
    output:
        "gs_marine_short.profile"
    shell:
        "wget https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/data/ground_truth/gs_marine_short.profile"

num_samples = int(config["num_samples"])
if num_samples < 1 or num_samples > 10:
    exit(-1)

rule combine_cami_results:
    input:
        expand("cami_result_sample{nums}.cami", nums=range(0, num_samples))
    output:
        "cami_result.cami"
    shell:
        "cat {input} > cami_result.cami"

rule run_cami:
    input:
        "cami_result.cami",
        "gs_marine_short.profile"
    output:
        "results/results.html"
    conda:
        "cami_env.yaml"
    shell:
        "opal.py -g {input[1]} {input[0]} -o ./results"
