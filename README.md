# ECS129-Seq-Comparision

Our implementation is based on Smale and colleagues' method for comparing protein sequences ["Towards a Mathematical Foundation of
Immunology and Amino Acid Chains"](https://arxiv.org/abs/1205.6031). We have created a console application in Python that can read in a BLOSUM62 text file, a beta parameter, and two sequence files in FASTA format. The program outputs K^3 hat and the distance between two sequences, which is a metric for the similarity between two sequences. Note that gaps are not considered as this program does not generate an alignment between two sequences. Additionally, it is difficult to compute a score between an amino acid and a gap (blank).

One way to run the `distance.py` program is as follows: 
```console
python3 distance.py blosum62matrix.txt 0.01 seq1.fa seq3.fa
```
