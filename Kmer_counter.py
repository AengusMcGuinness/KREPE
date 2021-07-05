#! /usr/bin/envpython3
import numpy as np
import os
import random
import math
import sys
import toyplot
import toyplot.browser
import toyplot.png
import matplotlib.pyplot as plt


kmer_length=int(sys.argv[1])
fasta_path=sys.argv[2]
plot_or_not=sys.argv[3]
de_bruijn_arg=sys.argv[4]
kmer_list=[]

def main():
    number_of_kmers = 0
    txt_path = fasta_path
    wrong_extension='.fasta'
    for character in wrong_extension:
        txt_path=txt_path.replace(character, "")
    base_path=txt_path
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
    if plot_or_not == 'plot_distribution':
        plotable_kmers = list(occurrence_dict.keys())
        frequency= values = list(occurrence_dict.values())
        plt.bar(range(len(occurrence_dict)), frequency, tick_label=plotable_kmers)
        plt.xlabel('Kmers')
        plt.ylabel('Frequency')
        plt.title('Distribution of Kmers')
        plt.show()
    elif plot_or_not == 'not_plot_distributiion':
        pass
    if de_bruijn_arg == 'plot_de_bruijn':
        edges = get_debruijn_edges_from_kmers(occurrence_dict)
        de_bruijn_graph = plot_debruijn_graph(edges)
        toyplot.png.render(de_bruijn_graph[0], f"debruijn_figure_{base_path}.png")
        toyplot.browser.show(de_bruijn_graph[0])
    elif de_bruijn_arg == 'not_plot_de_bruijn':
        pass

    
def get_debruijn_edges_from_kmers(kmers):
    """
    Every possible (k-1)mer (n-1 suffix and prefix of kmers) is assigned
    to a node, and we connect one node to another if the (k-1)mer overlaps 
    another. Nodes are (k-1)mers, edges are kmers.
    """
    # store edges as tuples in a set
    edges = set()
    
    # compare each (k-1)mer
    for k1 in kmers:
        for k2 in kmers:
            if k1 != k2:            
                # if they overlap then add to edges
                if k1[1:] == k2[:-1]:
                    edges.add((k1[:-1], k2[:-1]))
                if k1[:-1] == k2[1:]:
                    edges.add((k2[:-1], k1[:-1]))

    return edges

def plot_debruijn_graph(edges, width=800, height=800):
    graph=toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width = width,
        height = height,
        tmarker=">",
        vsize=25,
        vstyle={"stroke":"black", "stroke-width":1, "fill": "none"},
        vlstyle={"font-size": "12px"},
        estyle={"stroke": "black", "stroke-width":1},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    return graph



if __name__ == '__main__':
    main()

