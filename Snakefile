rule download_yacht:
    output:
	"env/yacht_env.yml"
    shell:
        "git clone git@github.com:KoslickiLab/YACHT.git && cd YACHT"

rule download_sample:
    output:
        "sample_file.tar.gz"
    shell:
        "wget -O {output} https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_0_reads.tar.gz"

rule unzip_sample_step1:
    input:
        "sample_file.tar.gz"
    output:
        "sample/simulation_short_read/2018.08.15_09.49.32_sample_0/reads/anonymous_reads.fq.gz",
        "sample/simulation_short_read/2018.08.15_09.49.32_sample_0/reads/reads_mapping.tsv.gz"
    shell:
        "mkdir -p sample && tar -xvzf {input} -C sample"

rule unzip_sample_step2:
    input:
        "sample/simulation_short_read/2018.08.15_09.49.32_sample_0/reads/anonymous_reads.fq.gz",
        "sample/simulation_short_read/2018.08.15_09.49.32_sample_0/reads/reads_mapping.tsv.gz"
    output:
        "sample/anonymous_reads.fq.gz",
        "sample/reads_mapping.tsv.gz"
    shell:
        "cp {input[0]} {output[0]} && cp {input[1]} {output[1]}"

rule unzip_sample_step3:
    input:
        "sample/anonymous_reads.fq.gz",
        "sample/reads_mapping.tsv.gz"
    output:
        "sample/anonymous_reads.fq",
        "sample/reads_mapping.tsv"
    shell:
        "gzip -dv {input[0]} && gzip -dv {input[1]}"

rule sketch_sample:
    input:
        "sample/anonymous_reads.fq"
    output:
        "sample.sig.zip"
    conda:
        "env/yacht_env.yml"
    shell:
        "yacht sketch sample --infile {input} --kmer 31 --scaled 3000 --outfile {output}"

rule unzip_ref:
    output:
        "refs/refdb_config.json"
    shell:
        "tar -xvzf gtdb_scale3000_pretrained.tar.gz"

rule run_yacht:
    input:
        "sample.sig.zip",
	"refs/refdb_config.json"
    output:
        "result.xlsx"
    conda:
        "env/yacht_env.yml"
    shell:
        "yacht run --json '{input[1]}' --sample_file '{input[0]}' --num_threads 32 --keep_raw --significance 0.99 --min_coverage_list 1 0.5 0.1 0.05 0.01 --out {output}"

rule create_genome_to_taxid:
    input:
        "result.xlsx"
    output:
        "genome_to_taxid.tsv"
    script:
        "create_genome_to_taxid.py"

rule convert_yacht_output:
    input:
        "result.xlsx",
        "genome_to_taxid.tsv"
    output:
        "cami_result.cami"
    conda:
        "env/yacht_env.yml"
    shell:
        "yacht convert --yacht_output '{input[0]}' --sheet_name 'raw_result' --genome_to_taxid '{input[1]}' --mode 'cami' --sample_name 'marmgCAMI2_short_read_sample_0' --outfile_prefix 'cami_result' --outdir ./"

rule download_gt:
    output:
        "gs_marine_short.profile"
    shell:
        "wget https://raw.githubusercontent.com/CAMI-challenge/second_challenge_evaluation/master/profiling/marine_dataset/data/ground_truth/gs_marine_short.profile"

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