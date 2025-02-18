import os
import sys
import glob
from Bio import SeqIO

def get_sequence_pair(fasta_file):
    """Extract reference and aligned sequences from a pairwise alignment file."""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    if len(sequences) != 2:
        print(f"Warning: {fasta_file} does not contain exactly 2 sequences")
        return None, None, None, None
        
    return str(sequences[0].seq), str(sequences[1].seq), sequences[0].id, sequences[1].id

def collect_sequence_pairs(fasta_files):
    """Collect all sequence pairs and find maximum length."""
    sequence_pairs = []
    max_length = 0
    ref_name = None
    
    for fasta_file in sorted(fasta_files):
        ref_seq, aln_seq, this_ref_name, aln_name = get_sequence_pair(fasta_file)
        if ref_seq is not None:
            if ref_name is None:
                ref_name = this_ref_name
            sequence_pairs.append((aln_name, ref_seq, aln_seq))
            max_length = max(max_length, len(ref_seq))
    
    return sequence_pairs, max_length, ref_name

def pad_sequence(seq, max_length):
    """Pad sequence with gaps to match max_length."""
    return seq + '-' * (max_length - len(seq))

def align_sequences(sequences):
    """
    Align sequences by adding gaps where needed to ensure all residues line up.
    Input: List of (sequence_name, reference_sequence, aligned_sequence) tuples
    Output: List of (sequence_name, aligned_reference, aligned_aligned) tuples
    """
    # Convert sequences to lists for easier manipulation
    seq_lists = [(name, list(ref), list(aln)) for name, ref, aln in sequences]
    max_length = max(len(ref) for _, ref, _ in sequences)
    
    # Process each position
    pos = 0
    while pos < max_length:
        # Check if any reference sequence has a gap at this position
        has_gap = any(pos < len(ref) and ref[pos] == '-' for _, ref, _ in seq_lists)
        
        if has_gap:
            # Insert gap in both reference and aligned sequences that have a residue at this position
            for _, ref, aln in seq_lists:
                if pos < len(ref) and ref[pos] != '-':
                    ref.insert(pos, '-')
                    aln.insert(pos, '-')
            max_length += 1
        pos += 1
    
    # Convert back to strings
    return [(name, ''.join(ref), ''.join(aln)) for name, ref, aln in seq_lists]

def process_alignments(sequence_pairs, max_length):
    """Return all sequences properly aligned."""
    padded_seqs = [(name, pad_sequence(ref, max_length), pad_sequence(aln, max_length)) 
                   for name, ref, aln in sequence_pairs]
    return align_sequences(padded_seqs)

def write_merged_fasta(aligned_sequences, output_file, ref_name):
    """Write aligned sequences to a merged FASTA file."""
    # Find maximum length across all sequences
    max_length = max(len(ref_seq) for _, ref_seq, _ in aligned_sequences)
    max_length = max(max_length, max(len(aln_seq) for _, _, aln_seq in aligned_sequences))
    
    with open(output_file, 'w') as f:
        # Write reference sequence once, padded to max length
        _, ref_seq, _ = aligned_sequences[0]
        padded_ref = pad_sequence(ref_seq, max_length)
        f.write(f">{ref_name}\n{padded_ref}\n")
        
        # Write all aligned sequences, padded to max length
        for seq_name, _, aln_seq in aligned_sequences:
            padded_aln = pad_sequence(aln_seq, max_length)
            f.write(f">{seq_name}\n{padded_aln}\n")

def extract_alignments(directory):
    """
    Extract and align sequences from all pairwise alignment files in a directory.
    
    Args:
        directory (str): Path to directory containing pairwise alignment FASTA files
    Returns:
        tuple: (aligned_sequences, reference_name) where aligned_sequences is a list of 
               tuples containing (sequence_name, reference_sequence, aligned_sequence)
    """
    fasta_files = glob.glob(os.path.join(directory, "*.fasta"))
    
    if not fasta_files:
        print(f"No FASTA files found in {directory}")
        return None, None
        
    sequence_pairs, max_length, ref_name = collect_sequence_pairs(fasta_files)
    return process_alignments(sequence_pairs, max_length), ref_name

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_directory> <output_file>")
        sys.exit(1)
        
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory")
        sys.exit(1)
        
    aligned_sequences, ref_name = extract_alignments(input_dir)
    
    try:
        os.remove(output_file)
    except:
        pass
    if aligned_sequences:
        write_merged_fasta(aligned_sequences, output_file, ref_name)

# if __name__ == '__main__':
#     main()
