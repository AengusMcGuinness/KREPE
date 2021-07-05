import random
import toyplot
import toyplot.browser
import toyplot.png

def main():
    binary = "0001110100"
    kmers= get_kmer_count_from_sequence(binary, k=3, cyclic=False)

    kmers
    binary= "0000110010111101"
    kmers= get_kmer_count_from_sequence(binary, k=4, cyclic=True)
    print(kmers)
    edges= get_debruijn_edges_from_kmers(kmers)
    edges

    print(binary)
    gr = plot_debruijn_graph(edges)
    toyplot.png.render(gr[0], "debuijn_figure.png")
    toyplot.browser.show(gr[0])
    


def get_kmer_count_from_sequence(sequence, k=3, cyclic=True):
    kmers = {}

    for i in range(0, len(sequence)):
        kmer = sequence[i:i +k]
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]
            else:
                if len(kmer) !=k:
                    continue
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
    return kmers


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


def plot_debruijn_graph(edges, width=500, height=500):
    graph=toyplot.graph(
        [i[0] for i in edges],
        [i[1] for i in edges],
        width = width,
        height = height,
        tmarker=">",
        vsize=25,
        vstyle={"stroke":"black", "stroke-width":2, "fill": "none"},
        vlstyle={"font-size": "11px"},
        estyle={"stroke": "black", "stroke-width":2},
        layout=toyplot.layout.FruchtermanReingold(edges=toyplot.layout.CurvedEdges()))
    return graph

main()
