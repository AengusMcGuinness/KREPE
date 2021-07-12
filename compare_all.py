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
import os
import random
import math
import sys
import matplotlib.pyplot as plt
from sourmash import fig, MinHash, SourmashSignature, save_signatures, load_one_signature
from itertools import combinations

kmer_length=int(sys.argv[1])
file_type = str(sys.argv[-1])
major_dict = {}
list_of_files = []
list_of_basenames = []
file_lengths = []
jaccards = []
sketch_list = []

def main():
    for k in sys.argv[2:-1]:
        major_dict.update({k : {}})
        formatted_file=file_cleaning(k)
        #print(formatted_file)
        character_count = len(formatted_file) - kmer_length
        file_lengths.append(character_count)
        #print(character_count)
        kmer_counting(character_count, formatted_file, k)
    byte_size = max(file_lengths)
    for k in sys.argv[2:-1]:
        sketch_list = implement_kmers(k, byte_size)
        print(sketch_list)
        jaccards_for_k = jaccardin_it_up(sketch_list)
        print(jaccards_for_k)        

def jaccardin_it_up(sketch_list):
    jaccards_for_k = []
    for a, b in combinations(sketch_list, 2):
        jaccard_index = a.jaccard(b)
        jaccards_for_k.append(jaccard_index)
    else:
        StopIteration
    jaccards.append(jaccards_for_k)
    return jaccards

def implement_kmers(k, byte_size):
    file_sketch = MinHash(n = byte_size, ksize = kmer_length)
    sketch_list.append(file_sketch)
    keys = major_dict[k].keys()
    for i in list(keys):
        incrementor = 0
        sketch_list[incrementor].add_kmer(i)
        incrementor += 1
    return sketch_list
        
def kmer_counting(character_count, formatted_file, k):
    #os.system('date --iso=seconds')
    for l in range(character_count):
        kmers= formatted_file[l:(l + kmer_length)]
        if not kmers in major_dict[k].keys():
            major_dict[k].update({kmers : 1})
        else:
             major_dict[k][kmers] += 1
        #os.system('date --iso=seconds')

def basename_creation(file_with_extension):
    if file_type == '-fasta' or '-fastq':
        basename = str(file_with_extension[:-6])
    if file_type == '-fna':
        basename = str(file_with_extension[:-4])
    list_of_basenames.append(basename)
    return basename

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
    if file_type == '-fasta' or '.fna':
        with open(unclean_file, 'r') as file:
            unformatted_file = file.readlines()
            for j in range(len(unformatted_file)):
                dna_lines = (j * 2) + 1
                if dna_lines > len(unformatted_file):
                    StopIteration
                else:
                     formatted_file=formatted_file  + unformatted_file[dna_lines]
    formatted_file = formatted_file.replace("\n", "")
    return formatted_file
    
main()


# def holy_list(unclean_file):
#     holy_list = []
#     os.system('date --iso=seconds')
#     formatted_file = ""
#     if file_type == '-fastq':
#         with open(unclean_file, 'r') as file:
#             unformatted_file = file.readlines()
#             for i in range(len(unformatted_file)):
#                 holy_number = (4*i + 1)
#                 if holy_number => len(unformatted_file):
#                     StopIteration
#                 else:
#                     holy_list.append(holy_number)
#         formatted_file=formatted_file+unformatted_file[for h in holy_list: h]                
# def file_cleaning(unclean_file):
#     os.system('date --iso=seconds')
#     formatted_file = ""
#     if file_type == '-fastq':
#         with open(unclean_file, 'r') as file:
#             unformatted_file = file.readlines()
#             for i in range(len(unformatted_file)):
#                 for j in range(len(unformatted_file)):
#                     if j == (1 + 4*i):
#                         formatted_file=unformatted_file[j] + formatted_file
#                     else:
#                         pass
#             formatted_file.replace("\n", "")
#             return formatted_file
#     else:
#         with open(unclean_file, 'r') as file:
#             unformatted_file = file.readlines()
#             for i in range(len(unformatted_file)):
#                 for j in range(len(unformatted_file)):
#                     if j == (i*2 + 1):
#                         formatted_file=unformatted_file[j] + formatted_file
#                     else:
#                         pass
#             formatted_file.replace("\n", "")
#             os.system('date --iso=seconds')
#             return formatted_file
#
# def if_removing_less_than_one():
#     for m in range(len(major_dict)):
#         if major_dict[k].values() < 2:
#             major_dict[k].pop()
# def jaccardin_it_up(sketches):
#     sketch_values=list(sketches.values())
#     incrementor = 0
#     for i in sketch_values:
#         incrementor += 1
#         if incrementor + 1 > len(sketch_values):
#             StopIteration
#         else:
#             jaccard = i.jaccard(sketch_values[incrementor + 1])
#             jaccards.append(jaccard)
