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

    # preprocess the BLOSUM 62 matrix
    # first row: order of amino acid along the column
    
    def process_bl62(bl62_file):  
        # open the file
        bl62_parse = open(bl62_file)
        next(bl62) # first line will just be the order of amino acids 

        # first letter on following rows, defines amino acid for that row
        # each letter contains a set of 20 scores  
        bl62 = [list(entry.split()) for entry in bl62_parse]
        bl62.close()
        
    # sequence 1 and sequence 2
    # https://stackoverflow.com/questions/18395587/splitting-characters-from-a-text-file-in-python
    def process_sequences(seq1_file, seq2_file):
        with open(seq1_file) as s1:
            seq1 = [list(line.rstrip()) for line in s1]

        with open(seq2_file) as s2:
            seq2 = [list(line.rstrip()) for line in s2]

        seq1.close()
        seq2.close()

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