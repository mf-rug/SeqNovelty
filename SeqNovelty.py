from Bio import AlignIO
from Bio.SeqIO import write
import json
import os
import re
import sys


def calculate_identity(ref_seq, query_seq, usable_columns, percent=True):
    """Calculate % identity between ref_seq and query_seq for specified columns."""
    matches = 0
    total = 0
    for idx in usable_columns:
        if ref_seq[idx] != "-" and query_seq[idx] != "-":
            total += 1
            if ref_seq[idx] == query_seq[idx]:
                matches += 1
    if percent:
        return (matches / total) * 100 if total > 0 else 0
    else:
        return matches

def find_best_hit(ref_seq, sequences, usable_columns):
    """Find the best matching sequence based on % identity."""
    best_hit_id = None
    best_identity = -1
    for seq_id, seq in sequences.items():
        identity = calculate_identity(ref_seq, seq, usable_columns, percent=False)
        if identity > best_identity:
            best_identity = identity
            best_hit_id = seq_id
    return best_hit_id, best_identity

def generate_alignment_visual(ref_seq, query_seq, matching_columns):
    """Generate a visual representation of the alignment with '|' for matches."""
    visual = ["|" if idx in matching_columns and ref_seq[idx] == query_seq[idx] else " " for idx in range(len(ref_seq))]
    return "".join(visual)

def iterative_alignment(input_fasta, ref_seq_id, output_dir="."):
    """Perform iterative alignment and best hit selection."""
    os.makedirs(output_dir, exist_ok=True)

    # Run initial alignment
    initial_alignment = os.path.join(output_dir, "initial_alignment.fasta")

    # Read alignment
    alignment = AlignIO.read(initial_alignment, "fasta")
    sequences = {record.description: str(record.seq) for record in alignment}
    ref_seq = sequences.pop(ref_seq_id)

    # Initialize variables
    usable_columns = set(range(len(ref_seq)))  # All columns initially usable
    results = {}

    iteration = 1
    used_sequences = []
    used_sequences_full = {}  # Store full sequences for all used hits
    summary = {
        ref_seq_id: {
            "matched_indices": []
        }
    }
    while sequences:
        # Find the best hit
        best_hit_id, best_identity = find_best_hit(ref_seq, sequences, usable_columns)
        if best_identity == 0 or best_hit_id is None:
            break

        best_hit_seq = sequences[best_hit_id]

        # Identify matching columns
        matching_columns = [
            idx for idx in usable_columns
            if ref_seq[idx] != "-" and best_hit_seq[idx] == ref_seq[idx]
            # if best_hit_seq[idx] == ref_seq[idx]
        ]

        # Remove matching columns from usable columns
        usable_columns -= set(matching_columns)
        summary[ref_seq_id]["matched_indices"].extend(matching_columns)
        
        # Store the full sequence and matched indices
        summary[best_hit_id] = {
            "matched_indices": matching_columns
        }
        used_sequences_full[best_hit_id] = best_hit_seq

        # Generate reference visual
        ref_visual = ["^" if idx not in usable_columns else " " for idx in range(len(ref_seq))]
        ref_visual = "".join(ref_visual)

        # Save pairwise alignment to a file
        pairwise_output = os.path.join(output_dir, f"pairwise_iter_{iteration}.fasta")
        with open(pairwise_output, "w") as f:
            f.write(f""">{ref_seq_id}
{ref_seq}
>{best_hit_id}
{best_hit_seq}
""")

        # Store results
        results[iteration] = {
            "hit_sequence": best_hit_id,
            "matching_columns": matching_columns,
            "percent_identity": best_identity,
        }

        # Track used sequences
        used_sequences.append((best_hit_id, matching_columns))
        used_sequences_full[best_hit_id] = best_hit_seq  # Keep full sequence for FASTA output

        # Remove best hit from sequences
        del sequences[best_hit_id]

        if len(sequences) == 0:
            print("All sequences have been analysed.")

        # Terminate if all residues are matched
        if not usable_columns:
            print("All reference residues have been matched.")
            break

        iteration += 1
    
    with open('summary.txt', "w") as f:
        ref_visual = re.sub(r'\^', '_', ref_visual)
        f.write(f">{ref_seq_id}\n{ref_visual}\n{ref_seq}")

        for seq_id, matching_columns in used_sequences:
            masked_sequence = [
                ref_seq[idx] if idx in matching_columns else " " for idx in range(len(ref_seq))
            ]
            f.write(f"\n>{seq_id}\n{''.join(masked_sequence)}")

    # Save the reference and matched sequences as a FASTA file
    fasta_output = os.path.join(output_dir, "final_alignment.fasta")
    with open(fasta_output, "w") as f:
        f.write(f">{ref_seq_id}\n{ref_seq}\n")
        for seq_id, _ in used_sequences:
            f.write(f">{seq_id}\n{used_sequences_full[seq_id]}\n")

    # Save results as a JSON file
    json_output = os.path.join(output_dir, "alignment_summary.json")
    with open(json_output, "w") as f:
        json.dump(summary, f, indent=4)

    # Save results to a JSON file
    results_file = os.path.join(output_dir, "results.json")
    with open(results_file, "w") as f:
        json.dump(results, f, indent=4)
    #print(f"Results saved to {results_file}.")


def main():
    input_fasta = sys.argv[1]  # Path to your input FASTA file
    seq_id = sys.argv[2]  # ID of the sequence to mask
    iterative_alignment(input_fasta, seq_id)
    print(f"SeqNovelty ran with {seq_id}, output saved to ./alignment_summary.json, ./results.json, ./final_alignment.fasta")


# if __name__ == '__main__':
#     main()