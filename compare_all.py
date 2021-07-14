'''

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
    os.system('date --iso=seconds')
    for k in sys.argv[2:-2]:
        labels.append(k)
        basename_creation(k)
        major_dict.update({k : {}})
        formatted_file=file_cleaning(k)
        character_count = len(formatted_file) - kmer_length
        file_lengths.append(character_count)
        kmer_counting(character_count, formatted_file, k)
        #system_stats()
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
            #labeltext.append(list_of_basenames[h]+entry)
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

def plot_composite_matrix(D, labeltext, show_labels=True, show_indices=True,
                          vmax=1.0, vmin=0.0, force=False):
    """Build a composite plot showing dendrogram + distance matrix/heatmap.
    Returns a matplotlib figure."""
    if D.max() > 1.0 or D.min() < 0.0:
        error('This matrix doesn\'t look like a distance matrix - min value {}, max value {}', D.min(), D.max())
        if not force:
            raise ValueError("not a distance matrix")
        else:
            notify('force is set; scaling to [0, 1]')
            D -= D.min()
            D /= D.max()

    if show_labels:
        show_indices = True

    fig = pylab.figure(figsize=(11, 8))
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])

    # plot dendrogram
    Y = sch.linkage(D, method='single')  # centroid

    dendrolabels = labeltext
    if not show_labels:
        dendrolabels = [str(i) for i in range(len(labeltext))]

    Z1 = sch.dendrogram(Y, orientation='left', labels=dendrolabels,
                        no_labels=not show_indices, get_leaves=True)
    ax1.set_xticks([])

    xstart = 0.45
    width = 0.45
    if not show_labels:
        xstart = 0.315
    scale_xstart = xstart + width + 0.01

    # re-order labels along rows, top to bottom
    idx1 = Z1['leaves']
    reordered_labels = [ labeltext[i] for i in reversed(idx1) ]

    # reorder D by the clustering in the dendrogram
    D = D[idx1, :]
    D = D[:, idx1]

    # show matrix
    axmatrix = fig.add_axes([xstart, 0.1, width, 0.6])

    im = axmatrix.matshow(D, aspect='auto', origin='lower',
                          cmap=pylab.cm.YlGnBu, vmin=vmin, vmax=vmax)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([scale_xstart, 0.1, 0.02, 0.6])
    pylab.colorbar(im, cax=axcolor)

    return fig, reordered_labels, D

def system_stats():
    print("Cpu Usage: ", psutil.cpu_percent())
    # gives an object with many fields
    print("Memory Usage: ", psutil.virtual_memory())
    # you can have the percentage of used RAM
    print("Memory Percentage: ", psutil.virtual_memory().percent)
    # you can calculate percentage of available memory
    print("Available Memory:", psutil.virtual_memory().available * 100 / psutil.virtual_memory().total)



main()

        #os.system('date --iso=seconds')
