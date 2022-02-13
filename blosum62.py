# ECS129 Project
# Emily Xiong, Betty Wu, and Wanting Zeng 
#
# Option 4: A Metric for Protein Sequences
#   Computing the similarity betweene two protein sequences
#   reads in: matrix BL62, parameter β, and 2 sequences (S1 and S2)
#   outputs the distance between two sequences


from os import linesep
import numpy as np  
import sys

# global variables
bl62 = {}
seq1 = {}
seq2 = {}
beta = 0 

# preprocess the BLOSUM 62 matrix
# matrix size 20 x 20
# Ref: https://www.biostars.org/p/405990/
def process_bl62(bl62_file):  
    bl62_file = open(bl62_file, 'r')
    lines = bl62_file.readlines()
    bl62_file.close()
    blosum = {} 
    aminoAcids = lines.pop(0) # first line: letters of the amino acids
    aminoAcids = aminoAcids.split() # get each letter
    
    for row in lines:
        entries = row.split()
        amino = entries.pop(0) 
        blosum[amino] = {}

        for col in aminoAcids:
            blosum[amino][col] = entries.pop(0)

    return(blosum)

# sequence 1 and sequence 2
# https://stackoverflow.com/questions/18395587/splitting-characters-from-a-text-file-in-python
# https://stackoverflow.com/questions/32953339/how-to-skip-line-with-matching-pattern-in-python?noredirect=1&lq=1
def process_sequences(sequence):
    seq_file = open(sequence, 'r')
    seq = {}
    seq = [list(line.rstrip()) for line in seq_file if not line.startswith(">")]
    seq_file.close()
    return seq
    # return [list(line.rstrip()) for line in seq if not line.startswith(">")]

# compute K1
    
# computes K2 given parameter K1

# computes K3 and distance given parameter K2

def main():
    # command line input
    # matrix BL62, parameter β, 2 sequences (s1, s2)
    if(len(sys.argv) != 5):
        sys.exit('Usage:[matrix BL62][parameter β][sequence 1][sequence 2]')
    
    bl62 = process_bl62(sys.argv[1])
    beta = sys.argv[2]
    seq1 = process_sequences(sys.argv[3])
    seq2 = process_sequences(sys.argv[4])
    # print(seq1)
    # print(seq2)
    # print(bl62)

# run the main program
main()