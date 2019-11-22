#!/usr/bin/env python3

import sys
import argparse
import resource

import numpy as np
import networkx as nx


import fasta.reader
import kmerator

def main():
  ap = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
  description="""
  guly.py, based on SWIGG.
  Given a list of fasta sequences, the script will produce a kmer graph in
  graphml chosen per various thresholds provided by user through arguments. This
  file can then be further visualized using graph visualizers such as Gephi or
  Cytoscope.""")

  ap.add_argument("-k", "--kmer-length",
                        type=int,
                        help="required, the length of k-mer that needs to be used",
                        required=True)
  ap.add_argument("-t", "--threshold",
                        type=int,
                        help="required, minimum threshold for kmers to be chosen\
                              \nshould be shared by at least these many sequences",
                        required=True)
  ap.add_argument("-f", "--fasta",
                        nargs="+",
                        help="required, set of all fasta sequences",
                        required=True)
  ap.add_argument("-rw", "--repeat_threshold_within",
                        type=int,
                        help="required, threshold for filtering kmers that can \
                              \nbe used for graph building within the same sequence. \
                              \nIf a certain kmer is present within the sequence \
                              \nabove the threshold, it will be discarded",
                        required=True)
  ap.add_argument("-ra", "--repeat_threshold_across",
                        type=int,
                        help="required, threshold for filtering kmers that can\
                              \nbe used for graph building across all the\
                              \nsequences we are looking at. If a certain kmer\
                              \nis chosen, but is present more than the threshold\
                              \nmentioned across all sequences, it will be discarded.\
                              \nThis is important to reduce repeats and avoid cyclic graphs.",
                      required=True)
  ap.add_argument("-o", "--out",
                      type=str,
                      default='kmer-graph.xml',
                      help="output path to prefix of output file name.",
                      required=True)
  ap.add_argument("-d", "--dumpkmer",
                         type=str,
                         default=None,
                         help="dump kmers in file for debugging. Can be large")
  args = ap.parse_args()
  ktor = kmerator.Kmerator(int(args.kmer_length), int(args.threshold), int(args.repeat_threshold_across))
  fp = fasta.reader.Reader()
  print("Searching {}-mers".format(ktor.kmer_len) , file=sys.stderr)
  for i in args.fasta:
    print(" in {}:".format(i), end = '', file=sys.stderr)
    fp.parse(i, ktor)
    print(" analyzed {} / {} kmers ( analyzed / distinct) in {} sequence(s)".format(ktor.kmer_count, ktor.db.kmer_count(),
                                                                                    ktor.db.sequence_count()), file=sys.stderr)
  #kmerator.filter(int(args.repeat_threshold_within))
  if args.dumpkmer:
    print("Dumping kmers into {}".format(args.dumpkmer), file=sys.stderr)
    ktor.dump_kmers(fout=args.dumpkmer)
  print("Assembling kmer graph into \'{}\'".format(args.out), file=sys.stderr)
  ktor.find_anchor_kmers(args.out)
  print("Ressource usage: {}".format(resource.getrusage(resource.RUSAGE_SELF)))
  return 0

if __name__ == '__main__':
  main()
