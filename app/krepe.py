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
import os
import math
import sys
import matplotlib.pyplot as plt
from sourmash import fig, MinHash
from matplotlib_venn import venn2
import re
import toyplot.pdf
import toyplot.browser
import toyplot
import matplotlib.pyplot as plt
from pylab import *

def venn_diagram(kmer_list_file_1, kmer_list_file_2, path1, path2):
    print("Number of Kmers File 1: ", len(kmer_list_file_1))
    print("Number of Kmers File 2: ", len(kmer_list_file_2))
    labels = [path1, path2]
    circle_one = set(kmer_list_file_1)
    circle_two = set(kmer_list_file_2)
    venn2([circle_one, circle_two], set_labels = labels)
    plt.show()
    
def genetic_distance_calc(jaccard_index, kmer_length):
    factor1 = -1/kmer_length
    factor2 = -1*(math.log(2 * jaccard_index) / 1 + jaccard_index)
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

def Kmers_of_file_1(formatted_file, kmer_length, occurrence_dict_file_1, kmer_list_file_1):
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
    key_list = list(occurrence_dict_file_1.keys())
    for key in key_list:
        if occurrence_dict_file_1[key] < 2:
            occurrence_dict_file_1.pop(key)
    return number_of_kmers_file_1

def Kmers_of_file_2(formatted_file, kmer_length, occurrence_dict_file_2,  kmer_list_file_2):
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
    key_list = list(occurrence_dict_file_2.keys())
    for key in key_list:
        if occurrence_dict_file_2[key] < 2:
            occurrence_dict_file_2.pop(key)
    return number_of_kmers_file_2
            
def get_debruijn_edges_from_kmers(kmers):
    edges = set()
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:            
                if k1[1:] == k2[:-1]:
                    edges.add((k1[:-1], k2[:-1]))
                if k1[:-1] == k2[1:]:
                    edges.add((k2[:-1], k1[:-1]))

    return edges

def basename_creation(file_with_extension, file_types, list_of_basenames):
    if file_types == '-fasta' or '-fastq':
        basename = str(file_with_extension[:-6])
    if file_types == '-fna':
        basename = str(file_with_extension[:-4])
    list_of_basenames.append(basename)
    return basename

def plot_debruijn_graph(edges, config_width, config_height, config_circle_width, config_line_width):
    graph=toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width = config_width,
        height = config_height,
        tmarker=">",
        vsize=25,
        vstyle={"stroke":"black", "stroke-width":config_circle_width, "fill": "none"},
        vlstyle={"font-size": "6px"},
        estyle={"stroke": "black", "stroke-width":config_line_width},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    return graph
        
def dendrogram(labeltext, sketch_list):
    jaccard_matrix = []
    jaccard_indexes = []
    for j in sketch_list:
        jaccard_indexes = []
        for l in sketch_list:
            index = j.jaccard(l)
            jaccard_indexes.append(index)
        jaccard_matrix.append(jaccard_indexes)
    matrix = np.matrix(jaccard_matrix)
    fig.plot_composite_matrix(matrix, labeltext)
    os.system('date --iso=seconds')
    plt.show()

def meta_data_aquisiton(meta_data, list_of_basenames):
    with open(meta_data, 'r') as f:
        meta_data_content = f.readlines()
        labeltext = []
        for h in range(len(list_of_basenames)):
            entry = meta_data_content[h].replace("\n", "")
            labeltext.append(entry)
    return labeltext
        
def implement_kmers(byte_size, k, major_dict, kmer_length, sketch_list):
    file_sketch = MinHash(byte_size, kmer_length)
    keys = major_dict[k].keys()
    for i in list(keys):
        file_sketch.add_kmer(i)
    sketch_list.append(file_sketch)
    return sketch_list
        
def kmer_counting(character_count, formatted_file, k, kmer_length, major_dict):
    for l in range(character_count):
        kmers=formatted_file[l:(l + kmer_length)]
        if not kmers in major_dict[k].keys():
            major_dict[k].update({kmers : 1})
        else:
             major_dict[k][kmers] += 1

def genome_visualization():
    kmer_length=int(sys.argv[1])
    path=sys.argv[2]
    plot_or_not=sys.argv[3]
    de_bruijn_arg=sys.argv[4]
    file_types=sys.argv[5]
    output = sys.argv[6]
    kmer_list=[]
    list_of_basenames = []
    number_of_kmers = 0
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
    occurrence_dict = {}
    character_count=len(formatted_file)-int(kmer_length)
    for i in range(character_count):
        number_of_kmers=number_of_kmers + 1
        kmers= formatted_file[i:(i + kmer_length)]
        if not kmers in occurrence_dict:
            occurrence_dict[kmers] = 1
        else:
            occurrence_dict[kmers] += 1
        kmer_list.append(kmers)
    if plot_or_not == 'plot_distribution':
        plotable_kmers = list(occurrence_dict.keys())
        frequency= values = list(occurrence_dict.values())
        plt.bar(range(len(occurrence_dict)), frequency, )
        plt.xlabel('Kmers')
        plt.ylabel('Frequency')
        plt.title('Distribution of Kmers')
        plt.show()
    elif plot_or_not == 'not_plot_distributiion':
        pass
    if de_bruijn_arg == 'plot_de_bruijn':
        de_bruijn_config = sys.argv[-1]
        with open(de_bruijn_config, 'r') as f:
            de_bruijn_parameters = f.readlines()
            config_width = int(de_bruijn_parameters[0].replace("\n", ""))
            config_height= int(de_bruijn_parameters[1].replace("\n", ""))
            config_circle_width= float(de_bruijn_parameters[2].replace("\n", ""))
            config_line_width= float(de_bruijn_parameters[3].replace("\n", ""))
        basename = basename_creation(path, file_types, list_of_basenames)
        edges = get_debruijn_edges_from_kmers(occurrence_dict)
        de_bruijn_graph = plot_debruijn_graph(edges, config_width, config_height, config_circle_width, config_line_width)
        toyplot.pdf.render(de_bruijn_graph[0], f"debruijn_figure_{basename}.pdf")
        toyplot.browser.show(de_bruijn_graph[0])
    elif de_bruijn_arg == 'not_plot_de_bruijn':
        pass
    if output == "-q":
        pass
    elif output == '-v':
        print(occurrence_dict)
    return

def organism_comparison():
    kmer_length=int(sys.argv[1])
    path1=sys.argv[2]
    path2=sys.argv[3]
    file_types=sys.argv[4]
    occurrence_dict_file_1 = {}
    occurrence_dict_file_2 = {}
    kmer_list_file_1 = []
    kmer_list_file_2 = []
    formatted_file = file_cleaning_two_point_o(path1)
    number_of_kmers_file_1=Kmers_of_file_1(formatted_file, kmer_length, occurrence_dict_file_1, kmer_list_file_1)
    formatted_file = file_cleaning_two_point_o(path2)
    number_of_kmers_file_2=Kmers_of_file_2(formatted_file, kmer_length, occurrence_dict_file_2,  kmer_list_file_2)
    if number_of_kmers_file_2 <= number_of_kmers_file_1:
        kmer_mash_file_2 = MinHash(len(kmer_list_file_1), ksize=kmer_length)    
        kmer_mash_file_1 = MinHash(len(kmer_list_file_1), ksize=kmer_length)
    else:
        kmer_mash_file_2 = MinHash(len(kmer_list_file_2), ksize=kmer_length)    
        kmer_mash_file_1 = MinHash(len(kmer_list_file_2), ksize=kmer_length)
    for i in kmer_list_file_1:
        kmer_mash_file_1.add_kmer(i)
    for i in kmer_list_file_2:
        kmer_mash_file_2.add_kmer(i)
    jaccard_index = kmer_mash_file_1.jaccard(kmer_mash_file_2)
    print("Jaccard Similarity: ", jaccard_index)
    genetic_distance = genetic_distance_calc(jaccard_index, kmer_length)
    print("Genetic Distance: ", genetic_distance)
    venn_diagram(kmer_list_file_1, kmer_list_file_2, path1, path2)
    return
    
def compare_all():
    kmer_length=int(sys.argv[1])
    meta_data=str(sys.argv[-1])
    file_types = str(sys.argv[-2])
    major_dict = {}
    labels = []
    list_of_basenames = []
    file_lengths = []
    sketch_list = []
    for k in sys.argv[2:-2]:
        labels.append(k)
        basename_creation(k, file_types, list_of_basenames)
        major_dict.update({k : {}})
        formatted_file = file_cleaning_two_point_o(k)
        character_count = len(formatted_file) - kmer_length
        file_lengths.append(character_count)
        kmer_counting(character_count, formatted_file, k, kmer_length, major_dict)
    byte_size = max(file_lengths)
    for k in sys.argv[2:-2]:
        sketch_list = implement_kmers(byte_size, k, major_dict, kmer_length, sketch_list)
    if meta_data == '-base':
        labeltext=labels
        dendrogram(labeltext, sketch_list)
    else:
        labeltext = meta_data_aquisiton(meta_data, list_of_basenames)
        dendrogram(labeltext, sketch_list)
    return
