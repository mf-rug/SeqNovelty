from Bio import AlignIO
import json
import os
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


def iterative_alignment_from_pairwise(alignment_dir, ref_seq_id):
    """Perform iterative alignment and best hit selection from pairwise alignments."""
    
    # Read all pairwise alignments and collect sequences
    sequences = {}
    ref_seq = None
    
    for filename in os.listdir(alignment_dir):
        if not filename.endswith('.fasta'):
            continue
        alignment = AlignIO.read(os.path.join(alignment_dir, filename), "fasta")
        for record in alignment:
            if record.description == ref_seq_id:
                if ref_seq is None:
                    ref_seq = str(record.seq).replace("-", "")  # Store ungapped version
                elif str(record.seq).replace("-", "") != ref_seq:
                    raise ValueError(f"Inconsistent reference sequence found in {filename}")
            else:
                sequences[record.description] = str(record.seq)

    if ref_seq is None:
        raise ValueError(f"Reference sequence {ref_seq_id} not found in alignments")

    # Initialize variables
    usable_columns = set(range(len(ref_seq)))
    results = {}
    used_sequences = []  # Will store (seq_id, sequence, matching_columns)
    summary = {ref_seq_id: {"matched_indices": []}}

    # Rest of the logic remains the same as iterative_alignment
    for iteration in range(1, len(sequences) + 1):
        best_hit_id, best_identity = find_best_hit(ref_seq, sequences, usable_columns)
        if best_identity == 0 or best_hit_id is None:
            break
        best_hit_seq = sequences.pop(best_hit_id)

        matching_columns = [
            idx for idx in usable_columns
            if ref_seq[idx] != "-" and best_hit_seq[idx] == ref_seq[idx]
        ]

        usable_columns -= set(matching_columns)
        summary[ref_seq_id]["matched_indices"].extend(matching_columns)
        summary[best_hit_id] = {"matched_indices": matching_columns}
        used_sequences.append((best_hit_id, best_hit_seq, matching_columns))

        results[iteration] = {
            "hit_sequence": best_hit_id,
            "matching_columns": matching_columns,
            "percent_identity": best_identity,
        }

        if not usable_columns:
            break

    return {
        "sequences": {
            "reference": {"id": ref_seq_id, "sequence": ref_seq},
            "hits": [(seq_id, sequence) for seq_id, sequence, _ in used_sequences]
        },
        "summary": summary,
        "results": results
    }

def iterative_alignment(initial_alignment, ref_seq_id):
    """Perform iterative alignment and best hit selection."""

    # Read alignment
    alignment = AlignIO.read(initial_alignment, "fasta")
    sequences = {record.description: str(record.seq) for record in alignment}
    ref_seq = sequences.pop(ref_seq_id)

    # Initialize variables
    usable_columns = set(range(len(ref_seq)))
    results = {}
    used_sequences = []  # Will store (seq_id, sequence, matching_columns)
    summary = {ref_seq_id: {"matched_indices": []}}

    for iteration in range(1, len(sequences) + 1):
        # Find the best hit
        best_hit_id, best_identity = find_best_hit(ref_seq, sequences, usable_columns)
        if best_identity == 0 or best_hit_id is None:
            break

        best_hit_seq = sequences.pop(best_hit_id)

        # Identify matching columns
        matching_columns = [
            idx for idx in usable_columns
            if ref_seq[idx] != "-" and best_hit_seq[idx] == ref_seq[idx]
        ]

        # Update data structures
        usable_columns -= set(matching_columns)
        summary[ref_seq_id]["matched_indices"].extend(matching_columns)
        summary[best_hit_id] = {"matched_indices": matching_columns}
        used_sequences.append((best_hit_id, best_hit_seq, matching_columns))

        # Store results
        results[iteration] = {
            "hit_sequence": best_hit_id,
            "matching_columns": matching_columns,
            "percent_identity": best_identity,
        }

        if not usable_columns:
            break

    # Remove file writing code and return the data instead
    return {
        "sequences": {
            "reference": {"id": ref_seq_id, "sequence": ref_seq},
            "hits": [(seq_id, sequence) for seq_id, sequence, _ in used_sequences]
        },
        "summary": summary,
        "results": results
    }

def main():
    input = sys.argv[1]  # Path to your input FASTA file
    seq_id = sys.argv[2]  # ID of the sequence to mask
    output_dir = "."
    
    if os.path.isfile(input):
        input_fasta = input
        alignment_data = iterative_alignment(input_fasta, seq_id)

    elif os.path.isdir(input):
        input_dir = input
        alignment_data = iterative_alignment_from_pairwise(input_dir, seq_id)
    
    else:
        print(f'Input {input} neither valid file nor directory')
        return

    # Get returned data from iterative_alignment
    # alignment_data = iterative_alignment(input_fasta, seq_id)
    
    # Write final alignment FASTA
    with open(os.path.join(output_dir, "final_alignment.fasta"), "w") as f:
        f.write(f">{alignment_data['sequences']['reference']['id']}\n{alignment_data['sequences']['reference']['sequence']}\n")
        for seq_id, sequence in alignment_data['sequences']['hits']:
            f.write(f">{seq_id}\n{sequence}\n")

    # Write JSON files
    for filename, data in [
        ("alignment_summary.json", alignment_data['summary']),
        ("results.json", alignment_data['results'])
    ]:
        with open(os.path.join(output_dir, filename), "w") as f:
            json.dump(data, f, indent=4)
            
    print(f"SeqNovelty ran with {seq_id}, output saved to ./alignment_summary.json, ./results.json, ./final_alignment.fasta")


# if __name__ == '__main__':
#     main()