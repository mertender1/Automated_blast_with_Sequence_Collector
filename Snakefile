import os
import subprocess

configfile: "config.yaml"

# --- PATHS ---
DB_FINAL = config["paths"]["db_fasta"]  
DB_DIR = os.path.dirname(DB_FINAL)
DB_FILTERED = os.path.join(DB_DIR, "filtered.faa")
DB_FILTERED_ADD = os.path.join(DB_DIR, "filtered_plus_extra.faa")
OUTDIR = config["paths"]["outdir"]

# --- GLOBAL VARS ---
MAX_ROUNDS = int(config["params"]["max_rounds"])
MAX_SEQS = int(config["params"]["max_seqs"])
QUERY = config["paths"]["initial_query"]
EXTRA_FAA = config["paths"].get("extra_faa")
HAS_EXTRA = EXTRA_FAA and os.path.exists(EXTRA_FAA)
FINAL_SOURCE = DB_FILTERED_ADD if HAS_EXTRA else DB_FILTERED

DB_TO_INDEX = DB_FINAL if config["params"].get("use_clustered_db", True) else FINAL_SOURCE

os.makedirs(DB_DIR, exist_ok=True)

# --- DYNAMIC HELPER FUNCTIONS ---
def count_fasta_seqs(fasta_file):
    if not os.path.exists(fasta_file): return 0
    try:
        res = subprocess.check_output(f"grep -c '>' {fasta_file}", shell=True)
        return int(res.decode().strip())
    except: return 0

def get_iterations(wildcards):
    """Evaluates checkpoints to decide if we need more rounds."""
    found_files = []
    for r in range(1, MAX_ROUNDS + 1):
        round_name = f"round_{r}"
        # This triggers the checkpoint for the current round
        checkpoint_output = checkpoints.mmseqs_cluster.get(round=round_name).output.fasta
        found_files.append(checkpoint_output)
        
        # Check if we hit the limit
        count = count_fasta_seqs(checkpoint_output)
        if count >= MAX_SEQS:
            return found_files
    return found_files

# --- RULES ---

rule all:
    input:
        DB_FINAL,
        f"{OUTDIR}/summary.tsv"

# 1. PRECISE SUBTRACTION & BATCH DOWNLOAD
rule download_and_filter_db:
    input:
        taxid_list = config["database"]["taxid_list_file"]
    output:
        full_fasta = os.path.join(DB_DIR, "full_unfiltered.faa"),
        filtered_fasta = DB_FILTERED
    params:
        excludes = config["database"]["exclude_taxids"],
        api_key = config["database"]["ncbi_api_key"],
        tmp_dir = os.path.join(DB_DIR, "tmp_download"),
        chunks_dir = os.path.join(DB_DIR, "chunks_fasta") # Ensure this is defined
    conda: "envs/blast_env.yaml"
    shell:
        """
        export NCBI_API_KEY={params.api_key}
        mkdir -p {params.tmp_dir}
        mkdir -p {params.chunks_dir}
        
        # 1. Build Query
        EX_QUERY="("
        for tid in {params.excludes}; do
            EX_QUERY="${{EX_QUERY}}txid${{tid}}[Organism:exp] OR "
        done
        EX_QUERY="${{EX_QUERY% OR }})"
        MAIN_TID=$(head -n 1 {input.taxid_list} | tr -d '[:space:]')
        
        # 2. Fetch Accessions (Skip if file exists to save time)
        if [ ! -s {params.tmp_dir}/clean_ids.txt ]; then
            echo "Querying NCBI for accessions..."
            esearch -db protein -query "txid${{MAIN_TID}}[Organism:exp] NOT ${{EX_QUERY}}" | \
            efetch -format acc > {params.tmp_dir}/clean_ids.txt
        fi
        
        TOTAL_IDS=$(wc -l < {params.tmp_dir}/clean_ids.txt)
        echo "Found ${{TOTAL_IDS}} clean accessions."

        # 3. Split the file into chunks of 500 lines
        mkdir -p {params.tmp_dir}/id_chunks
        split -l 500 {params.tmp_dir}/clean_ids.txt {params.tmp_dir}/id_chunks/batch_
        
        # 4. Process each chunk file with header count validation
        BATCHES=$(ls {params.tmp_dir}/id_chunks/batch_*)
        for batch_file in $BATCHES; do
            batch_name=$(basename "$batch_file")
            id_list=$(tr '\\n' ',' < "$batch_file" | sed 's/,$//')
            EXPECTED_COUNT=$(wc -l < "$batch_file")
            
            success=false
            attempt=1
            while [ "$success" = false ] && [ $attempt -le 3 ]; do
                # Download (Overwrites if previous attempt was partial)
                efetch -db protein -id "$id_list" -format fasta > {params.chunks_dir}/seq_${{batch_name}}.fasta
                
                # Count headers (">") to verify integrity
                ACTUAL_COUNT=$(grep -c ">" {params.chunks_dir}/seq_${{batch_name}}.fasta || echo 0)
                
                if [ "$ACTUAL_COUNT" -ge "$EXPECTED_COUNT" ]; then
                    success=true
                else
                    echo "Batch ${{batch_name}} mismatch: Got ${{ACTUAL_COUNT}}, expected ${{EXPECTED_COUNT}}. Retry $attempt/3..."
                    sleep 5
                    attempt=$((attempt+1))
                fi
            done

            if [ "$success" = false ]; then
                echo "CRITICAL: Batch ${{batch_name}} failed after 3 attempts. Exiting."
                exit 1
            fi
            
            sleep 0.3
        done

        # 5. All batches verified. Combine them now.
        echo "--- Combining chunks into {output.full_fasta} ---"
        cat {params.chunks_dir}/seq_*.fasta > {output.full_fasta}

        echo "--- Deduplicating ---"
        seqkit rmdup -s {output.full_fasta} -o {output.filtered_fasta}
        
        # Cleanup
        rm -rf {params.tmp_dir}/id_chunks
        rm -rf {params.chunks_dir}
        """

# 1.5 MERGE EXTRA SEQUENCES
# This rule only runs if EXTRA_FAA is actually defined and exists
rule merge_extra_faa:
    input:
        filtered = DB_FILTERED,
        extra = EXTRA_FAA if EXTRA_FAA else []
    output:
        combined = DB_FILTERED_ADD
    shell:
        "cat {input.filtered} {input.extra} > {output.combined}"

# 2. INITIAL CLUSTERING (CD-HIT)
rule cluster_db_initial:
    input:
        # If EXTRA_FAA is provided, use the merged file. Otherwise, use the download.
        FINAL_SOURCE
    output: 
        DB_FINAL
    params:
        c = config["params"]["cdhit_id"]
    threads: 4
    conda: "envs/blast_env.yaml"
    shell: "cd-hit -i {input} -o {output} -c {params.c} -n 5 -T {threads} -M 0"

# 3. INDEXING
rule index_db:
    input: fasta = DB_TO_INDEX
    output: multiext(DB_TO_INDEX, ".pin", ".phr", ".psq")
    conda: "envs/blast_env.yaml"
    shell: "makeblastdb -in {input.fasta} -dbtype prot -parse_seqids"

# 4. BLAST ROUNDS
rule run_blast:
    input:
        query = lambda wildcards: QUERY if wildcards.round == "round_1" 
                else f"{OUTDIR}/round_{int(wildcards.round.split('_')[1])-1}/clustered.faa",
        db_files = multiext(DB_TO_INDEX, ".pin", ".phr", ".psq")
    output:
        tabular = f"{OUTDIR}/{{round}}/blast_results.tab"
    params:
        db = DB_TO_INDEX,
        evalue = config["params"]["blast"]["evalue"],
        outfmt = config["params"]["blast"]["outfmt"],
        threads = config["params"]["blast"]["threads"]
    conda: "envs/blast_env.yaml"
    shell:
        """
        blastp -query {input.query} -db {params.db} -out {output.tabular} \
               -evalue {params.evalue} -outfmt "{params.outfmt}" \
               -num_threads {params.threads}
        """

# 5. EXTRACTION
rule extract_hits:
    input: 
        tabular = f"{OUTDIR}/{{round}}/blast_results.tab",
        db_files = multiext(DB_TO_INDEX, ".pin", ".phr", ".psq") 
    output: 
        fasta = f"{OUTDIR}/{{round}}/extracted_hits.faa"
    params: 
        db = DB_TO_INDEX
    conda: "envs/blast_env.yaml"
    script: "scripts/run_blast.py"

# 6. MMSEQS CLUSTERING
checkpoint mmseqs_cluster:
    input: fasta = f"{OUTDIR}/{{round}}/extracted_hits.faa"
    output: fasta = f"{OUTDIR}/{{round}}/clustered.faa"
    params:
        tmp = f"{OUTDIR}/{{round}}/tmp_mmseqs",
        prefix = f"{OUTDIR}/{{round}}/mmseqs_out",
        min_id = config["params"]["mmseqs"]["identity"],
        cov = config["params"]["mmseqs"]["coverage"],
        sens = config["params"]["mmseqs"]["sensitivity"],
        cov_mode = config["params"]["mmseqs"]["cov_mode"]
    conda: "envs/blast_env.yaml"
    threads: 4
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs easy-cluster {input.fasta} {params.prefix} {params.tmp} \
            --min-seq-id {params.min_id} -c {params.cov} --cov-mode {params.cov_mode} \
            --cluster-mode 2 -s {params.sens} --threads {threads}
        
        [ -f {params.prefix}_rep_seq.fasta ] && mv {params.prefix}_rep_seq.fasta {output.fasta} || \
        ([ -f {params.prefix}_rep_seqs.fasta ] && mv {params.prefix}_rep_seqs.fasta {output.fasta} || cp {input.fasta} {output.fasta})
        """

# 7. SUMMARY
rule generate_summary:
    input:
        # Link to the checkpoint evaluator
        clustered_list = get_iterations
    output: 
        f"{OUTDIR}/summary.tsv"
    run:
        with open(output[0], "w") as f:
            f.write("Round\tRaw_Hits\tClustered_Reps\n")
            for clu_path in input.clustered_list:
                round_dir = os.path.dirname(clu_path)
                raw_path = os.path.join(round_dir, "extracted_hits.faa")
                
                raw_cnt = count_fasta_seqs(raw_path)
                clu_cnt = count_fasta_seqs(clu_path)
                
                r_label = os.path.basename(round_dir)
                f.write(f"{r_label}\t{raw_cnt}\t{clu_cnt}\n")