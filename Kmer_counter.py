#! /usr/bin/envpython3
import numpy as np
import os
import random
import math
import sys
import matplotlib.pyplot as plt


kmer_length=int(sys.argv[1])
fasta_path=sys.argv[2]
plot_or_not=sys.argv[3]
kmer_list=[]

def main():
    number_of_kmers = 0
    txt_path = fasta_path
    wrong_extension='.fasta'
    for character in wrong_extension:
        txt_path=txt_path.replace(character, "")
    cmd = f"grep -v 'length=' {fasta_path} > {txt_path}.txt" 
    os.system(cmd)
    print(f'{txt_path}')
    occurrence_dict = {}
    with open(f'{txt_path}.txt', 'r') as file:
        os.system('date --iso=seconds')
        txt_file = file.read().replace('\n', '')
        os.system('date --iso=seconds')
        character_count=len(txt_file)-int(kmer_length)
        for i in range(character_count):
            number_of_kmers=number_of_kmers + 1
            kmers= txt_file[i:(i + kmer_length)]
            if not kmers in occurrence_dict:
                occurrence_dict[kmers] = 1
            else:
                occurrence_dict[kmers] += 1
            kmer_list.append(kmers)
        os.system('date --iso=seconds')
    if len(occurrence_dict.keys()) < 100:
        print(occurrence_dict)
    else:
        print('[very big occurrence dict]')
    print(4**kmer_length, len(occurrence_dict.keys()))
    # list_of_keys_with_more_than_one_occurrence = []
    key_list = list(occurrence_dict.keys())
    for key in key_list:
        if occurrence_dict[key] < 2:
            # list_of_keys_with_more_than_one_occurrence.append(key)
            occurrence_dict.pop(key)
    print('======================================================')
    if len(occurrence_dict.keys()) < 100:
        print(occurrence_dict)
    else:
        print('[very big occurrence dict]')
    print(4**kmer_length, len(occurrence_dict.keys()))
    if plot_or_not == 'plot':
        plotable_kmers = list(occurrence_dict.keys())
        frequency= values = list(occurrence_dict.values())
        plt.bar(range(len(occurrence_dict)), frequency, tick_label=plotable_kmers)
        plt.xlabel('Kmers')
        plt.ylabel('Frequency')
        plt.title('Distribution of Kmers')
        plt.show()
    elif plot_or_not == 'not_plot':
        pass
        
if __name__ == '__main__':
    main()

