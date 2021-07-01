from Bio.Seq import Seq
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import math
import pandas as pd
import struct
import mmh3 
import bitarray
import sys
#These are packages I think we might need so I am just putting them at the top

'''-----------------------------Work in Progress---------------------------------------------'''

#Working on bloom filter class
random.seed(1)
def main():
    incrementer={}
    kmer_list=[]
    number_of_kmers=0
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
            temporary_data_structure.update(data_entry)
        file.close()
    fp = float(input("Enter a false probability as a decimal: "))
    print(temporary_data_structure)



def bloomfilter():
   # inserted_kmers=len(set(kmer_list))
    total_bits=math.ceil(number_of_kmers*(1.44*(math.log(fp, 2))))
    bloomeyfilter = bitarray(total_bits)
    bloomeyfilter.setall(0)
    number_of_hash=(total_bits / number_of_kmers) * np.log(2)
    #Not yet applicable
    #bloomfilter = BloomFilter(number_of_kmers, fp)


def hashing(kmer_list, number_of_hash):
    for i in (kmer_list):
       # mmh3.hash(kmer_list[i], random.randint(0, number_of_hash))
        hash_key=int(mmh3.hash(kmer_list[i]))
        f'{hash_key}incrementer' = 0
        if bloomeyfilter[hash_key::(hash_key + 1)] == bitarray('0'):
            del bloomeyfilter[hash_key::(hash_key + 1)]
            bitarray('1') = bloomeyfilter[hash_key::(hash_key + 1)]
        else:
            f'{hash_key}incrementer' += 1 


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
