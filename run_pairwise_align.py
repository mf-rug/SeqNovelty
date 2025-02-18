from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
import os
import sys
import glob

def calculate_identity(seq1, seq2):
    """Calculate the percentage identity between two aligned sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    total = len(seq1)
    return (matches / total) * 100

def create_pairwise_alignments(input_fasta, ref_seq_id, identity_threshold):
    """
    Create pairwise alignments between a reference sequence and all other sequences in a FASTA file.
    
    Args:
        input_fasta (str): Path to input FASTA file containing multiple sequences
        ref_seq_id (str): ID of the reference sequence to align against
        identity_threshold (float): Minimum percent identity to write the alignment
    """
    # Create output directory if it doesn't exist
    output_dir = "./pairwise_alignments"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # remove old pairwise files
    fasta_files = glob.glob(os.path.join(output_dir, "*.fasta"))
    for file in fasta_files:
        os.remove(file)

    # Read all sequences from input FASTA
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Find reference sequence
    ref_seq = None
    other_seqs = []
    for seq in sequences:
        if seq.id == ref_seq_id:
            ref_seq = seq
        else:
            other_seqs.append(seq)
            
    if ref_seq is None:
        raise ValueError(f"Reference sequence {ref_seq_id} not found in {input_fasta}")
        
    # Create aligner object with BLOSUM62 matrix and free end gaps
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10  # Standard gap opening penalty
    aligner.extend_gap_score = -0.5  # Standard gap extension penalty
    aligner.target_end_gap_score = 0  # Free end gaps
    aligner.query_end_gap_score = 0  # Free end gaps

    # Create pairwise alignments
    for seq in other_seqs:
        # Perform pairwise alignment
        alignments = aligner.align(ref_seq.seq, seq.seq)
        
        # Get best alignment
        best_alignment = alignments[0]
        
        aligned_ref = best_alignment[0]
        aligned_query = best_alignment[1]
        
        # Calculate percent identity
        identity = calculate_identity(aligned_ref, aligned_query)
        
        # Skip writing if identity is below the threshold
        if identity < identity_threshold:
            print(f"Skipping {seq.id} due to low identity: {identity:.2f}%")
            continue
        else:
            print(f"Accepting {seq.id} due to {identity:.2f}% > {identity_threshold}")
        
        # Write alignment to output file
        output_file = os.path.join(output_dir, f"{ref_seq_id}_vs_{seq.id}.fasta")
        with open(output_file, "w") as f:
            f.write(f">{ref_seq.id}\n{aligned_ref}\n")
            f.write(f">{seq.id}\n{aligned_query}\n")

def main():
    input_fasta = sys.argv[1]
    ref_seq_id = sys.argv[2]
    identity_threshold = float(sys.argv[3])  # Get the identity threshold from command line
    create_pairwise_alignments(input_fasta=input_fasta, ref_seq_id=ref_seq_id, identity_threshold=identity_threshold)
    print(f"run_pairwise ran with {ref_seq_id} on {input_fasta}, outputs saved to ./pairwise_alignments/")



# if __name__ == '__main__':
#     main()