# ECS129 Project
# Emily Xiong, Betty Wu, and Wanting Zeng 
#
# Option 4: A Metric for Protein Sequences
#   Computing the similarity betweene two protein sequences
#   reads in: matrix BL62, parameter β, and 2 sequences (S1 and S2)
#   outputs the distance between two sequences

# Last Update: 2/12/22 11:22 PM

import numpy as np  
import sys

# global variables
bl62 = {}
seq1 = {}
seq2 = {}
beta = 0.0

# computes K2 given parameter K1

# computes K3 and distance given parameter K2


# preprocess the BLOSUM 62 matrix and gets the K1 value 
# Ref: https://www.biostars.org/p/405990/
def process_bl62(bl62_file, beta):  
    bl62_file = open(bl62_file, 'r')
    lines = bl62_file.readlines()
    bl62_file.close()
    blosum62 = {} 

    # first line: letters of the amino acids
    aminoAcids = lines.pop(0)   
    aminoAcids = aminoAcids.split()     
    
    # separate each amino acid
    for row in lines:   # first letter defines amino acid for the row
        entries = row.split()   
        amino = entries.pop(0)  
        blosum62[amino] = {}

        # get score for each amino acid
        for col in aminoAcids:  # BL62(x,y)
            entry = entries.pop(0)
            blosum62[amino][col] = pow(float(entry), beta) # K1(x,y) = BL62(x,y)^β

    return np.array(blosum62)   # return blosum62

# sequence 1 and sequence 2
# https://stackoverflow.com/questions/18395587/splitting-characters-from-a-text-file-in-python
# https://stackoverflow.com/questions/32953339/how-to-skip-line-with-matching-pattern-in-python?noredirect=1&lq=1
# https://stackoverflow.com/questions/4978787/how-to-split-a-string-into-a-list-of-characters
# https://www.w3schools.com/python/ref_func_range.asp
def process_sequences(seq_file):
    seq_file = open(seq_file, 'r')
    lines = seq_file.readlines()
    seq_file.close()
    letters = ""
    
    for i in range(0, len(lines)):
        lines[i] = lines[i].strip('\n')
        line = lines[i]
        if line[i] == '>':  # skip 
            continue
        else:
            letters += line
    sequence = list(letters)
    return sequence 
    # return np.array(sequence)
    # return [list(line.rstrip()) for line in seq if not line.startswith(">")]
    
def main():
    # command line input
    # matrix BL62, parameter β, 2 sequences (s1, s2)
    if(len(sys.argv) != 5):
        sys.exit('Usage:[matrix BL62][parameter β][sequence 1][sequence 2]')
    
    beta = float(sys.argv[2])
    bl62 = process_bl62(sys.argv[1], beta)
    seq1 = process_sequences(sys.argv[3])
    seq2 = process_sequences(sys.argv[4])

    # used for testing
    print("Sequence1: ", seq1, "\n")
    print("Sequence2: ", seq2, "\n")
    print("Blosum62: ", bl62, "\n")

# run the main program
main()