# ECS129 Project
# Betty Wu, Emily Xiong, and Wanting Zeng
#
# Option 4: A Metric for Protein Sequences
#   Computing the similarity between two protein sequences
#   reads in: matrix BL62, parameter β, and 2 sequences (S1 and S2)
#   outputs the distance between two sequences
# 
# Special Thanks to Professor Koehl for spotting bugs

import math
import sys
import os.path
import random
from time import process_time
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

# global variables for sharing between functions
k1 = {}
aminoAcids = []

# @compute_k2 and @compute_k3 are versions that didn't utilize memoization
def compute_k2(u,v,k):
    product = 1
    for i in range(k):
        product*=k1[u[i]][v[i]]
    return product

def compute_k3(f, g):
    sum = 0
    min_len = min(len(f), len(g))
    # k is seq length
    for k in range(1,min_len+1):
        for i in range(len(f)-k+1):
            for j in range(len(g)-k+1):
                sum+=compute_k2(f[i:i+k], g[j:j+k],k)

    return sum

def speedup_compute_k3(f, g):

    sum = 0
    k2_memo={}

    min_len = min(len(f), len(g))
    for k in range(1,min_len+1):
        for i in range(len(f)-k+1):
            for j in range(len(g)-k+1):
                u = f[i:i+k]
                v = g[j:j+k]
                last_k1 = k1[u[k-1]][v[k-1]]
                if k == 1:
                    k2_memo[(k,i,j)] = last_k1
                else:
                    k2_memo[(k,i,j)]=k2_memo[(k-1,i,j)]*last_k1
                sum+= k2_memo[(k,i,j)]

    return sum

def normalized_k3(f,g):
    normalized = speedup_compute_k3(f,g)/math.sqrt(speedup_compute_k3(f,f)*speedup_compute_k3(g,g))
    print(f"Normalized K3: {normalized:.4f}")
    return normalized

def compute_distance(f, g):
    return math.sqrt(2*(1-normalized_k3(f, g)))

# preprocess the BLOSUM 62 matrix and gets the K1 value
def process_bl62(bl62_file, beta):
    if os.path.exists(bl62_file) is False:
         sys.exit("Error: BLOSUM 62 matrix file does not exist")
    bl62_file = open(bl62_file, 'r')
    lines = bl62_file.readlines()
    bl62_file.close()

    # first line: letters of the amino acids
    global aminoAcids
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
def process_sequences(seq_filename):
    if os.path.exists(seq_filename) is False:
         sys.exit("Error: sequence file " + seq_filename + " does not exist")

    seq_fd = open(seq_filename, 'r')
    lines = seq_fd.readlines()
    seq_fd.close()
    letters = ""

    for i in range(0, len(lines)):
        lines[i] = lines[i].strip('\n')
        line = lines[i]
        if line[i] == '>':  # skip
            continue
        else:
            letters += line
    sequence = list(letters.upper())

    # check if the sequence is valid
    global aminoAcids
    for amino in sequence:
        if amino not in aminoAcids:
            sys.exit(f"Error: sequence file {seq_filename} contains invalid amino acid {amino}")

    return sequence

def plot_shuffled_distr(k3hat, seq1, seq2, n):
    print("Start shuffling!")
    # our computed k3hat would be the first element
    results = [k3hat]
    for i in range(n):
        print(f"{i} round")
        shuffled_seq1 = random.sample(seq1, len(seq1))
        shuffled_seq2 = random.sample(seq2, len(seq2))
        shuffled_k3hat = normalized_k3(shuffled_seq1, shuffled_seq2)
        results.append(shuffled_k3hat)
    sns.kdeplot(results)
    data = np.array(results)
    zs = stats.zscore(data)
    print(f"The normalized k3 {data[0]:.4f} has z score {zs[0]:.4f}")
    plt.show()

def plot_runtime_for_fragments(seq1, seq2, step):
    results=[]
    max_fragment_size = min(len(seq1), len(seq2))
    if step > max_fragment_size:
        sys.exit(f"Error: step size {step} is too large. Try using different step size.")

    possible_indices = list(np.arange(step,max_fragment_size))
    end_indices = possible_indices[::step] # end index
    if max_fragment_size-1 not in end_indices:
        end_indices.append(max_fragment_size-1)
    for i in end_indices:
        timer_start = process_time()
        compute_distance(seq1[:i], seq2[:i])
        timer_end = process_time()
        results.append(timer_end-timer_start)
    sns.regplot(x=end_indices, y=results)
    plt.show()

def main():
    # command line input
    # BL62, parameter β, 2 sequences (s1, s2)
    if(len(sys.argv) != 5):
        sys.exit('Usage:[matrix BL62][parameter β][sequence 1][sequence 2]')    
    
    global k1

    beta = float(sys.argv[2])
    if(beta < 0):
        sys.exit('Error: beta must be positive') 
    k1 = process_bl62(sys.argv[1], beta)
    seq1 = process_sequences(sys.argv[3])
    seq2 = process_sequences(sys.argv[4])

    timer_start = process_time()
    distance = compute_distance(seq1, seq2)
    timer_end = process_time()

    print(f"Distance between two sequences: {distance:.4f}")
    print(f"Time elapsed: {timer_end-timer_start}")
    # plot_shuffled_distr(normalized_k3(seq1, seq2), seq1, seq2, 1000)
    # plot_runtime_for_fragments(seq1, seq2, 5)

# run the main program
if __name__ == '__main__':
    main()


