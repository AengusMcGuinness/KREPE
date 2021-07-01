from Bio.Seq import Seq
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import math
import pandas as pd
import struct
import mmh3 
from bitarray import bitarray
import sys
#These are packages I think we might need so I am just putting them at the top

'''-----------------------------Work in Progress---------------------------------------------'''

#Working on bloom filter class
random.seed(1)
def main():
    kmer_list=[]
    number_of_kmers = 0
    temporary_data_structure={}
    kmer_length=int(input("Enter Kmer Lengths: "))
    fasta_path=input("File path (FASTA FORMAT): ")
    txt_path = fasta_path
    wrong_extension='.fasta'
    for character in wrong_extension:
        txt_path=txt_path.replace(character, "")
    cmd = f"grep -v '[^[:alpha:]]' {fasta_path} > {txt_path}.txt" 
    os.system(cmd)
    with open(f'{txt_path}.txt', 'r') as file:
        txt_file = file.read().replace('\n', '')
        #print(txt_file)
        character_count=len(txt_file)-int(kmer_length)
        for i in range(character_count):
            number_of_kmers=number_of_kmers + 1
            kmers= txt_file[i:(i + kmer_length)]
            #data_entry={str(hash(kmers)):str(kmers)}
            kmer_list.append(kmers)
   #         temporary_data_structure.update(data_entry)
        file.close()
    fp = float(input("Enter a false probability as a decimal: "))
  #  print(temporary_data_structure)
    hash_funcs, bloomeyfilter = bloomfilter(number_of_kmers, fp)
    #print(hash_funcs)
    hashing(kmer_list, hash_funcs, bloomeyfilter, total_bits)
    #print(kmer_counting_dictionary)
    #print(number_of_kmers, total_bits, kmer_counting_dictionary)
    #print(total_bits)
    print(kmer_counting_dictionary)
    
def bloomfilter(number_of_kmers, fp):
   # inserted_kmers=len(set(kmer_list))
    global total_bits
    total_bits=int(abs(math.ceil(number_of_kmers*(1.44*(math.log(fp, 2))))))
    print(total_bits)
    bloomeyfilter = bitarray(total_bits)
    bloomeyfilter.setall(0)
    hash_funcs=math.ceil(total_bits / number_of_kmers) * np.log(2)
    return hash_funcs, bloomeyfilter
    #Not yet applicable
    #bloomfilter = BloomFilter(number_of_kmers, fp)


def hashing(kmer_list, number_of_hash, bloomeyfilter, total_bits):
    global kmer_counting_dictionary
    kmer_counting_dictionary={}
   # print(kmer_list)
    for i in range(len(kmer_list)):
       # mmh3.hash(kmer_list[i], random.randint(0, number_of_hash))
      #  hash_key=int((abs(mmh3.hash(kmer_list[(i)])))/10000)
        random_hash=int(random.randint(0, total_bits))
     #   print(total_bits)
     #   print(hash_key)
       # print(f'{kmer_list[i]}incrementer')
        kmer_counting_dictionary.update({kmer_list[i]: 0})
        if (bloomeyfilter[random_hash:(random_hash + 1)]) == bitarray('0'):
           bloomeyfilter[random_hash:(random_hash + 1)] = bitarray('1')
           # fill_one = bloomeyfilter[hash_key::(hash_key + 1)]
        else:
           kmer_counting_dictionary[kmer_list[i]] += 1
    print(bloomeyfilter)
    return kmer_counting_dictionary

# #txt_file= str(open(f'{txt_path}.txt', 'w'))
#     #now the counting begins
#     #blomfilter = BloomFilter()


# '''
if __name__ == '__main__':
    main()

        
# Useful for dealing with bits and filiping them, uses struct packages, probablt needs to be tinkered on, those functions at the bottom are from a previous project I have done

# def float_to_bin(num):
#     """Given a float, return a string with the individual bits"""
#     result = format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')
#     return result

# def bin_to_float(binary):
#     return struct.unpack('!f',struct.pack('!I', int(binary, 2)))[0]
#python already has a hash function so we can use that with our bloom filter

# def bitflip(x, pos):
#      fs = pack('f',x)
#      bval = list(unpack('BBBB',fs))
#      [q,r] = divmod(pos,8)
#      bval[q] ^= 1 << r
#      fs = pack('BBBB', *bval)
#      fnew=unpack('f',fs)
#      return fnew[0]
