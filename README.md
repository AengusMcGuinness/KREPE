# KREPE

## BACKGROUND INFORMATION:
In bioinformatics, it's important to be able to use strands of DNA
to compare organisms, understand and assemble genomes, and get a better
picture of how everything works together in an organism. A common tool
is the k-mer. When a biologist is working on figuring out a genome, or
comparing a cat and a dog,they start with the k-mers. The way they work
is similar to a paper shredder; you put in a paper (or multiple reads of
the same sequence), and shred it into bits. Then, you take all the paper
bits and piece them back together in the same way as before, by finding
similarities in each of the pieces. This is useful to biologists because
sometimes there are errors when the machine reads a sequence, so with
multiple copies it's simple to find errors and toss them out early so they
don't affect the overall 'picture'.

These k-mers are sometimes strung together to make stringy looking
graphs called 'De Bruijn Graphs', that are used to find a pathway through
a genome- this process of finding a pathway is called collapsing. 'Collapsing'
a De Bruijn graph is comparable to cobbling together different bus or train
routes that start at the place where the last one ends, and finding a way
to your destination.



## DESCRIPTION:
 The goal of the project is to create a streamlined, easy to use
 package and command line tool that allows users to count k-mers
 (nucleotide sequences of length k) and easily generates the
 following:

 - A list of the occurring k-mers in the given file that occur more
   than once
 - A histogram that shows the frequency of said k-mer occurrence
 - A De Bruijn graph that can then be used for genome assembly

 The package is protected by copyright under the GNU General Public
 License

## INSTALLATION:
 Use the package manager [pip]
 (https://pip.pypa.io/en/stable/) to install krepe

 ```bash
 pip install krepe
 ```
OR you can install without pip

```bash
$ git clone https://github.com/RGBwalnut/Kmer-Counting-Analysis
```

## USAGE:

 The inputs that are allowed are .fasta, .fastq, and .fna. 

## ROADMAP:
 This package is completely open-sourced and is being updated by
 authors Erika Pedersen and Aengus McGuinness as needs arise

## CONTRIBUTING:
 Contributions will be accepted, in the instance that the contributor
 is willing to conform with the following:

 A. be listed as a minor contributor in the description

 B. not be listed under the copyright

 C. be willing to hear suggestions on minor changes to the proposed
 contribution As well as in the instance the authors deem the
 contribution useful.


