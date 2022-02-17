# ECS129 Project
# Emily Xiong, Betty Wu, and Wanting Zeng
#
# Option 4: A Metric for Protein Sequences
#   Computing the similarity between two protein sequences
#   reads in: matrix BL62, parameter β, and 2 sequences (S1 and S2)
#   outputs the distance between two sequences

# Last Update: 2/17/22 3:27PM

from cmath import sqrt
import numpy as np
import sys

MAX_K = 3
k1 = {}

# computes K2 given parameter K1
def compute_k2(u,v,k):
    product = 1

    for i in range(k):
        #print("ui is ",u[i])
        #print("vi is ", v[i])
        #print("score is ", k1[u[i]][v[i]])
        product*=k1[u[i]][v[i]]
    return product

# computes K3 and distance given parameter K2
def compute_k3(f, g):
    #print("f is ", f)
    #print("g is ", g)

    sum = 0
    min_len = min(len(f), len(g))
    # k is seq length
    for k in range(1,MAX_K+1):
        #print("k is ",k)
        for i in range(min_len-k+1):
            #print("i is ",i)
            #print("f sub is ", f[i:i+k])
            #print("g sub is ", g[i:i+k])
            #print()
            sum+=compute_k2(f[i:i+k], g[i:i+k],k)
    return sum

def normalized_k3(f,g):
    return compute_k3(f,g)/sqrt(compute_k3(f,f)*compute_k3(g,g))

def compute_distance(f, g):
    return sqrt(2*(1-normalized_k3(f, g)))

# preprocess the BLOSUM 62 matrix and gets the K1 value
# Ref: https://www.biostars.org/p/405990/
def process_bl62(bl62_file, beta):
    bl62_file = open(bl62_file, 'r')
    lines = bl62_file.readlines()
    bl62_file.close()

    # first line: letters of the amino acids
    aminoAcids = lines.pop(0)
    aminoAcids = aminoAcids.split()

    # separate each amino acid
    for row in lines:   # first letter defines amino acid for the row
        entries = row.split()
        amino = entries.pop(0)
        k1[amino] = {}

        # get score for each amino acid
        for col in aminoAcids:  # BL62(x,y)
            bl62_entry = entries.pop(0)
            k1[amino][col] = pow(
                float(bl62_entry), beta)  # K1(x,y) = BL62(x,y)^β
    return k1

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


def main():
    # command line input
    # matrix BL62, parameter β, 2 sequences (s1, s2)
    if(len(sys.argv) != 5):
        sys.exit('Usage:[matrix BL62][parameter β][sequence 1][sequence 2]')

    global k1
    global MAX_K

    beta = float(sys.argv[2])
    k1 = process_bl62(sys.argv[1], beta)
    seq1 = process_sequences(sys.argv[3])
    seq2 = process_sequences(sys.argv[4])
    MAX_K = min(len(seq1), len(seq2))

    distance = compute_distance(seq1, seq2)
    #distance = compute_distance(seq1[:3], seq2[:3])


    # used for testing
    #print("Sequence1: ", seq1, "\n")
    #print("Sequence2: ", seq2, "\n")
    #print("Blosum62: ", k1, "\n")
    print("computed distance: ",distance)


# run the main program
if __name__ == '__main__':
    main()
