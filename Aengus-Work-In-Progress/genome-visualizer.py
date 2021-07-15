#! /usr/bin/envpython3
import numpy as np
import os
import random
import math
import sys
import toyplot
import toyplot.browser
import toyplot.png
import toyplot.pdf
import matplotlib.pyplot as plt
import re

kmer_length=int(sys.argv[1])
path=sys.argv[2]
plot_or_not=sys.argv[3]
de_bruijn_arg=sys.argv[4]
file_types=sys.argv[5]
output = sys.argv[6]
kmer_list=[]

def main():
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
        basename = basename_creation(path)
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

def basename_creation(file_with_extension):
    if file_types == '-fasta' or '-fastq':
        basename = str(file_with_extension[:-6])
    if file_types == '-fna':
        basename = str(file_with_extension[:-4])
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

if __name__ == '__main__':
    main()
