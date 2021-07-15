#! /usr/bin/envpython3
import numpy as np
import os
import random
import math
import sys
import matplotlib.pyplot as plt
from sourmash import MinHash
from matplotlib_venn import venn2
import re

kmer_length=int(sys.argv[1])
path1=sys.argv[2]
path2=sys.argv[3]
file_types=sys.argv[4]
occurrence_dict_file_1 = {}
occurrence_dict_file_2 = {}
kmer_list_file_1 = []
kmer_list_file_2 = []

def main():
    formatted_file = file_cleaning_two_point_o(path1)
    number_of_kmers_file_1=Kmers_of_file_1(formatted_file)
    formatted_file = file_cleaning_two_point_o(path2)
    number_of_kmers_file_2=Kmers_of_file_2(formatted_file)
    if number_of_kmers_file_2 <= number_of_kmers_file_1:
        print('hashing')
        kmer_mash_file_2 = MinHash(len(kmer_list_file_1), ksize=kmer_length)    
        kmer_mash_file_1 = MinHash(len(kmer_list_file_1), ksize=kmer_length)
    else:
        print('hashing')
        kmer_mash_file_2 = MinHash(len(kmer_list_file_2), ksize=kmer_length)    
        kmer_mash_file_1 = MinHash(len(kmer_list_file_2), ksize=kmer_length)
    for i in kmer_list_file_1:
        kmer_mash_file_1.add_kmer(i)
        print(i)
    for i in kmer_list_file_2:
        kmer_mash_file_2.add_kmer(i)
        print(i)
    jaccard_index = kmer_mash_file_1.jaccard(kmer_mash_file_2)
    print("Jaccard Similarity: ", jaccard_index)
    genetic_distance = genetic_distance_calc(jaccard_index, kmer_length)
    print("Genetic Distance: ", genetic_distance)
    venn_diagram()


def venn_diagram():
    #set_labels=(f'{path1}', f'{path2}')
    print()
    print(len(kmer_list_file_1))
    print(len(kmer_list_file_2))
    print(kmer_list_file_1)
    print(kmer_list_file_2)
    labels = [path1, path2]
    circle_one = set(kmer_list_file_1)
    circle_two = set(kmer_list_file_2)
    venn2([circle_one, circle_two], set_labels = labels)
    plt.show()
    
def genetic_distance_calc(jaccard_index, kmer_length):
    factor1 = -1/kmer_length
    factor2 = math.log(2 * jaccard_index) / 1 + jaccard_index
    genetic_distance = factor1 * factor2
    return genetic_distance

def file_cleaning_two_point_o(path):
    formatted_file = []
    nucleotides = re.compile('[ACTG]')
    with open(path, 'r') as f:
        unformatted_file = f.readlines()
        for j in range(len(unformatted_file)):
            if len(''.join((nucleotides.findall(unformatted_file[j])))) == (len((unformatted_file[j]).replace("\n", ""))):
                formatted_file.append(unformatted_file[j].replace("\n", ""))
            else:
                pass
    formatted_file = ''.join(formatted_file)
    return formatted_file


def Kmers_of_file_1(formatted_file):
    number_of_kmers_file_1= 0
    character_count_file_1=len(formatted_file)-int(kmer_length)
    for i in range(character_count_file_1):
        number_of_kmers_file_1 += 1
        kmers= formatted_file[i:(i + kmer_length)]
        if not kmers in occurrence_dict_file_1:
            occurrence_dict_file_1[kmers] = 1
        else:
            occurrence_dict_file_1[kmers] += 1
            kmer_list_file_1.append(kmers)
    os.system('date --iso=seconds')
    key_list = list(occurrence_dict_file_1.keys())
    for key in key_list:
        if occurrence_dict_file_1[key] < 2:
            occurrence_dict_file_1.pop(key)
    return number_of_kmers_file_1


def Kmers_of_file_2(formatted_file):
    number_of_kmers_file_2 = 0
    character_count_file_2=len(formatted_file)-int(kmer_length)
    for i in range(character_count_file_2):
        number_of_kmers_file_2 += 1
        kmers= formatted_file[i:(i + kmer_length)]
        if not kmers in occurrence_dict_file_2:
            occurrence_dict_file_2[kmers] = 1
        else:
            occurrence_dict_file_2[kmers] += 1
            kmer_list_file_2.append(kmers)
    os.system('date --iso=seconds')
    key_list = list(occurrence_dict_file_2.keys())
    for key in key_list:
        if occurrence_dict_file_2[key] < 2:
            occurrence_dict_file_2.pop(key)
    return number_of_kmers_file_2


    
def print_count(occurrence_dict):
    print('======================================================')
    if len(occurrence_dict.keys()) < 100:
        print(occurrence_dict)
    else:
        print('[very big occurrence dict]')
        print(4**kmer_length, len(occurrence_dict.keys()))



# def file_creation(path):
#     txt_path = path
#     if file_types == '-fastq':
#         txt_path = txt_path.replace(txt_path[-6:], "")
#         cmd_fastq = f"sed -n '1~4s/^@/>/p;2~4p' {path} | grep -v 'length=' > {txt_path}.txt"
#         os.system(cmd_fastq)
#     if file_types == '-fna':
#         txt_path = txt_path.replace(txt_path[-4:], "")
#         cmd = f"grep -v 'length=' {path} | grep -v '>' > {txt_path}.txt" 
#         os.system(cmd)
#     elif file_types == '-fasta':
#         txt_path = txt_path.replace(txt_path[-6:], "")
#         cmd = f"grep -v 'length=' {path} | grep -v '>' > {txt_path}.txt" 
#         os.system(cmd)
#     return txt_path
        
if __name__ == '__main__':
    main()
