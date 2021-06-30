from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
import struct

#These are packages I think we might need so I am just putting them at the top

#-----------------------------Work in Progress-----------------------------------------------

#I don't think the bloom filter needs to be pandas dataframe so I am just leaving it as a dictionary for now.
bloom_filter={} 

def main():
    





#python already has a hash function so we can use that with our bloom filter
def hashing():
    




    
if __name__ == "__main__":
    main()

# Useful for dealing with bits and filiping them, uses struct packages, probablt needs to be tinkered on, those functions at the bottom are from a previous project I have done

# def bitflip(x, pos):
#     fs = pack('f',x)
#     bval = list(unpack('BBBB',fs))
#     [q,r] = divmod(pos,8)
#     bval[q] ^= 1 << r
#     fs = pack('BBBB', *bval)
#     fnew=unpack('f',fs)
#     return fnew[0]

# def float_to_bin(num):
#     """Given a float, return a string with the individual bits"""
#     result = format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')
#     return result

# def bin_to_float(binary):
#     return struct.unpack('!f',struct.pack('!I', int(binary, 2)))[0]
