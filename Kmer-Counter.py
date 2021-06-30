from Bio.Seq import Seq
import numpy as np
import os
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

def main():
    data = []
    kmer_length=input("Enter Kmer Lengths: ")
    fasta_path=input("File path (FASTA FORMAT): ")
    txt_path = file_path
    wrong_extension='.fasta'
    for character in wrong_extension:
        txt_path=txt_path.replace(character, "")
    cmd = f"grep -v '[^[:alpha:]]' {file_path} > {txt_path}.txt" 
    os.system(cmd)
    txt_file= open(f'{txt_path}.txt', 'r')
    
        
    #now the counting begins
    #blomfilter = BloomFilter()

# Currently trying to understand this. Found this from a tutorial

# class BloomFilter():

#     '''
#     Class for Bloom filter, using murmur3 hash function
#     '''
#     def __init__(self, items_count, fp_prob):
#         '''
#         items_count : int
#             Number of items expected to be stored in bloom filter
#         fp_prob : float
#             False Positive probability in decimal
#         '''
#         # False possible probability in decimal
#         self.fp_prob = fp_prob

#         # Size of bit array to use
#         self.size = self.get_size(items_count, fp_prob)

#         # number of hash functions to use
#         self.hash_count = self.get_hash_count(self.size, items_count)

#         # Bit array of given size
#         self.bit_array = bitarray(self.size)

#         # initialize all bits as 0
#         self.bit_array.setall(0)
        

#     def add(self, item):
#         '''
#         Add an item in the filter
#         '''
#         digests = []
#         for i in range(self.hash_count):

#             # create digest for given item.
#             # i work as seed to mmh3.hash() function
#             # With different seed, digest created is different
#             digest = mmh3.hash(item, i) % self.size
#             digests.append(digest)

#             # set the bit True in bit_array
#             self.bit_array[digest] = True

#     def check(self, item):
#         '''
#         Check for existence of an item in filter
#         '''
#         for i in range(self.hash_count):
#             digest = mmh3.hash(item, i) % self.size
#             if self.bit_array[digest] == False:

#                 # if any of bit is False then,its not present
#                 # in filter
#                 # else there is probability that it exist
#                 return False
#         return True

#     @classmethod
#     def get_size(self, n, p):
#         '''
#         Return the size of bit array(m) to used using
#         following formula
#         m = -(n * lg(p)) / (lg(2)^2)
#         n : int
#             number of items expected to be stored in filter
#         p : float
#             False Positive probability in decimal
#         '''
#         m = -(n * math.log(p))/(math.log(2)**2)
#         return int(m)

#     @classmethod
#     def get_hash_count(self, m, n):
#         '''
#         Return the hash function(k) to be used using
#         following formula
#         k = (m/n) * lg(2)

#         m : int
#             size of bit array
#         n : int
#             number of items expected to be stored in filter
#         '''
#         k = (m/n) * math.log(2)
#         return int(k)




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
