'''
   This program counts the kmers in 2 files returns a kmer sketch as well as genetic distance betweent the input files.
   Copyright (C) 2021 Aengus McGuinness and Erika Pedersen
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
'''
import numpy as np
import pylab
import os
import sys
import math
from sourmash import MinHash, fig
import psutil

kmer_length=int(sys.argv[1])
file_type = str(sys.argv[-1])
major_dict = {}
labels = []
sketch_list = []
list_of_basenames = []
file_lengths = []

def main():
    os.system('date --iso=seconds')
    for k in sys.argv[2:-1]:
        labels.append(k)
        major_dict.update({k : {}})
        formatted_file=file_cleaning(k)
        character_count = len(formatted_file) - kmer_length
        file_lengths.append(character_count)
        kmer_counting(character_count, formatted_file, k)
    byte_size = max(file_lengths)
    for k in sys.argv[2:-1]:
        sketch_list = implement_kmers(byte_size, k)

def dendrogram():
    jaccard_matrix = []
    jaccard_indexes = []
    for j in sketch_list:
        jaccard_indexes = []
        for l in sketch_list:
            jaccard_indexes.append(j.jaccard(l))
        print(jaccard_indexes)
        jaccard_matrix.append(jaccard_indexes)
    matrix = np.matrix(jaccard_matrix)
    fig.plot_composite_matrix(matrix, labels)
    os.system('date --iso=seconds')
    plt.show()

def implement_kmers(byte_size, k):
    file_sketch = MinHash(byte_size, kmer_length)
    keys = major_dict[k].keys()
    for i in list(keys):
        file_sketch.add_kmer(i)
    sketch_list.append(file_sketch)
    return sketch_list

def kmer_counting(character_count, formatted_file, k):
    for l in range(character_count):
        kmers= formatted_file[l:(l + kmer_length)]
        if not kmers in major_dict[k].keys():
            major_dict[k].update({kmers : 1})
        else:
             major_dict[k][kmers] += 1

def file_cleaning(unclean_file):
    formatted_file = ""
    if file_type == '-fastq':
        with open(unclean_file, 'r') as file:
            unformatted_file = file.readlines()
            for j in range(len(unformatted_file)):
                dna_lines = (j * 4) +1
                if dna_lines > len(unformatted_file):
                    StopIteration
                else:
                    formatted_file=formatted_file  + unformatted_file[dna_lines]
    # if file_type == '-fasta' or '.fna':
    #     with open(unclean_file, 'r') as file:
    #         unformatted_file = file.readlines()
    #         for j in range(len(unformatted_file)):
    #             dna_lines = (j * 2) + 1
    #             if dna_lines > len(unformatted_file):
    #                 StopIteration
    #             else:
    #                  formatted_file=formatted_file  + unformatted_file[dna_lines]
    # formatted_file = formatted_file.replace("\n", "")
    return formatted_file

main()
