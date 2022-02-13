# ECS129 Project
# Emily Xiong, Betty Wu, and Wanting Zeng 
#
# Option 4: A Metric for Protein Sequences
#   Computing the similarity betweene two protein sequences
#   reads in: matrix BL62, parameter β, and 2 sequences (S1 and S2)
#   outputs the distance between two sequences

from matplotlib.pyplot import close
import numpy as np  
import sys

class Blosum62:

    # preprocess the input and store appropiately (matrix)
    def process(bl62_file, beta_par, seq1_file, seq2_file):  
        bl62 = np.loadtxt(bl62_file, dtype=int)
        beta = beta_par
        seq1 = np.loadtxt(seq1_file)
        seq2 = np.loadtxt(seq2.file)
        bl62.close()
        seq1.close()
        seq2.close()
    # Ref: https://stackoverflow.com/questions/29629315/import-text-file-as-matrix-in-numpy

    # compute K1
        

# computes K2 given parameter K1

# computes K3 and distance given parameter K2

def main():
    # command line input
    # matrix BL62, parameter β, 2 sequences (s1, s2)
    if(len(sys.argv) != 5):
        sys.exit('Usage:[matrix BL62][parameter β][sequence 1][sequence 2]')

    bl62_file = sys.argv[1]
    beta_par = sys.argv[2]
    seq1_file = sys.argv[3]
    seq2_file = sys.argv[4]
    blosum62 = Blosum62(bl62_file, beta_par, seq1_file, seq2_file)

# run the main program
main()