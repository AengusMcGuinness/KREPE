#! /usr/bin/envpython3
import numpy as np
import os
import random
import math
import sys
import matplotlib.pyplot as plt
from sourmash import MinHash

kmer_length=int(sys.argv[1])
fasta_path1=sys.argv[2]
fasta_path2=sys.argv[3]
occurrence_dict_file_1 = {}
occurrence_dict_file_2 = {}
kmer_list_file_1 = []
kmer_list_file_2 = []

def main():
    txt_path = file_creation(fasta_path1)
    number_of_kmers_file_1=Kmers_of_file_1(f"{txt_path}.txt")
    txt_path = file_creation(fasta_path2)
    number_of_kmers_file_2=Kmers_of_file_2(f"{txt_path}.txt")
    if number_of_kmers_file_2 <= number_of_kmers_file_1:
        kmer_mash_file_2 = MinHash(number_of_kmers_file_1, ksize=kmer_length)     
        kmer_mash_file_1 = MinHash(number_of_kmers_file_1, kmer_length)
    else:
        kmer_mash_file_2 = MinHash(number_of_kmers_file_2, ksize=kmer_length)     
        kmer_mash_file_1 = MinHash(number_of_kmers_file_2, kmer_length)
    for i in kmer_list_file_1:
        kmer_mash_file_1.add_kmer(i)
    for i in kmer_list_file_2:
        kmer_mash_file_2.add_kmer(i)
    print(kmer_mash_file_1.jaccard(kmer_mash_file_2))

def file_creation(path):
    txt_path = path
    txt_path = txt_path.replace(txt_path[-6:], "")
    cmd = f"grep -v 'length=' {path} > {txt_path}.txt" 
    os.system(cmd)
    return txt_path

def Kmers_of_file_1(txt_path):
    with open(txt_path, 'r') as file:
        os.system('date --iso=seconds')
        number_of_kmers_file_1= 0
        txt_file = file.read().replace('\n', '')
        os.system('date --iso=seconds')
        character_count_file_1=len(txt_file)-int(kmer_length)
        for i in range(character_count_file_1):
            number_of_kmers_file_1 += 1
            kmers= txt_file[i:(i + kmer_length)]
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


def Kmers_of_file_2(txt_path):
    with open(txt_path, 'r') as file:
        os.system('date --iso=seconds')
        txt_file = file.read().replace('\n', '')
        os.system('date --iso=seconds')
        character_count_file_2=len(txt_file)-int(kmer_length)
        number_of_kmers_file_2 = 0
        for i in range(character_count_file_2):
            number_of_kmers_file_2 += 1
            kmers= txt_file[i:(i + kmer_length)]
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
    if len(occurrence_dict.keys()) < 100:
        print(occurrence_dict)
    else:
        print('[very big occurrence dict]')
        print(4**kmer_length, len(occurrence_dict.keys()))
    print('======================================================')
    if len(occurrence_dict.keys()) < 100:
        print(occurrence_dict)
    else:
        print('[very big occurrence dict]')
        print(4**kmer_length, len(occurrence_dict.keys()))
            
if __name__ == '__main__':
    main()
