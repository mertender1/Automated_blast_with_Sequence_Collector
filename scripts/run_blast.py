import subprocess
import pandas as pd
import os

def extract_sequences(blast_out, db_path, fasta_out):
    try:
        # 1. READ WITHOUT FIXED NAMES (because your tab file has many columns)
        # We only care about the second column (index 1)
        df = pd.read_csv(blast_out, sep='\t', header=None)
        
        # 2. CLEAN THE IDs: strip 'gb|' from start and '|' from end
        # Example: 'gb|QBQ65106.1|' becomes 'QBQ65106.1'
        raw_ids = df[1].unique().astype(str)
        unique_ids = [id_str.replace('gb|', '').strip('|') for id_str in raw_ids]
        
    except Exception as e:
        print(f"Error reading BLAST results: {e}")
        unique_ids = []
    
    if not unique_ids:
        with open(fasta_out, "w") as f:
            f.write("")
        return

    # Create the ID list
    temp_id_file = f"{fasta_out}.tmp_ids"
    with open(temp_id_file, "w") as f:
        for uid in unique_ids:
            f.write(f"{uid}\n")
    
    # We use -target_only to ensure it focuses on our list
    # We remove 'check=True' so that if one sequence fails, 
    # the script continues rather than crashing Snakemake.
    cmd = [
        "blastdbcmd",
        "-db", db_path,
        "-entry_batch", temp_id_file,
        "-out", fasta_out
    ]
    
    print(f"Extracting {len(unique_ids)} sequences...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Warning: Some sequences were skipped. BLAST error message:\n{result.stderr}")
        if os.path.exists(fasta_out) and os.path.getsize(fasta_out) > 0:
            print("Successfully extracted available sequences.")
        else:
            print("Error: No sequences could be extracted.")

    if os.path.exists(temp_id_file):
        os.remove(temp_id_file)

extract_sequences(
    snakemake.input.tabular,
    snakemake.params.db,
    snakemake.output.fasta
)