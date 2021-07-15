#! /usr/bin/envpython3
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
from pylab import *
import os
import math
import sys
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from sourmash import fig, MinHash, SourmashSignature
from matplotlib_venn import venn3
import psutil


kmer_length=int(sys.argv[1])
meta_data=str(sys.argv[-1])
file_type = str(sys.argv[-2])
major_dict = {}
labels = []
sketch_list = []
list_of_basenames = []
file_lengths = []

def main():
    for k in sys.argv[2:-2]:
        labels.append(k)
        basename_creation(k)
        major_dict.update({k : {}})
        formatted_file = file_cleaning_two_point_o(k)
        character_count = len(formatted_file) - kmer_length
        file_lengths.append(character_count)
        kmer_counting(character_count, formatted_file, k)
    byte_size = max(file_lengths)
    for k in sys.argv[2:-2]:
        sketch_list = implement_kmers(byte_size, k)
    if meta_data == '-base':
        labeltext=labels
        dendrogram(labeltext)
    else:
        labeltext = meta_data_aquisiton()
        dendrogram(labeltext)
        
def dendrogram(labeltext):
    jaccard_matrix = []
    jaccard_indexes = []
    for j in sketch_list:
        jaccard_indexes = []
        for l in sketch_list:
            jaccard_indexes.append(j.jaccard(l))
        jaccard_matrix.append(jaccard_indexes)
    matrix = np.matrix(jaccard_matrix)
    fig.plot_composite_matrix(matrix, labeltext)
    os.system('date --iso=seconds')
    plt.show()

def meta_data_aquisiton():
    with open(meta_data, 'r') as f:
        meta_data_content = f.readlines()
        labeltext = []
        for h in range(len(list_of_basenames)):
            entry = meta_data_content[h].replace("\n", "")
            labeltext.append(entry)
    return labeltext

def basename_creation(file_with_extension):
    if file_type == '-fasta' or '-fastq':
        basename = str(file_with_extension[:-6])
    if file_type == '-fna':
        basename = str(file_with_extension[:-4])
    list_of_basenames.append(basename)
    return basename
        
def implement_kmers(byte_size, k):
    file_sketch = MinHash(byte_size, kmer_length)
    keys = major_dict[k].keys()
    for i in list(keys):
        file_sketch.add_kmer(i)
    sketch_list.append(file_sketch)
    return sketch_list
        
def kmer_counting(character_count, formatted_file, k):
    for l in range(character_count):
        kmers=formatted_file[l:(l + kmer_length)]
        if not kmers in major_dict[k].keys():
            major_dict[k].update({kmers : 1})
        else:
             major_dict[k][kmers] += 1

def file_cleaning_two_point_o(unclean_file):
    formatted_file = []
    nucleotides = re.compile('[ACTG]')
    with open(unclean_file, 'r') as f:
        unformatted_file = f.readlines()
        for j in range(len(unformatted_file)):
            if len(''.join((nucleotides.findall(unformatted_file[j])))) == (len((unformatted_file[j]).replace("\n", ""))):
                formatted_file.append(unformatted_file[j].replace("\n", ""))
            else:
                pass
    formatted_file = ''.join(formatted_file)
            #print(formatted_file)
    return formatted_file    

def system_stats():
    print("Cpu Usage: ", psutil.cpu_percent())
    # gives an object with many fields
    print("Memory Usage: ", psutil.virtual_memory())
    # you can have the percentage of used RAM
    print("Memory Percentage: ", psutil.virtual_memory().percent)
    # you can calculate percentage of available memory
    print("Available Memory:", psutil.virtual_memory().available * 100 / psutil.virtual_memory().total)

if __name__ == '__main__':
    main()
