#!/usr/bin/env python3
#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------

from . import sequence

class Reader:

  @staticmethod
  def parse(fname, kmerator):
    fh = open(fname, 'r')
    header = None
    seq = ''
    for i in fh:
      if i[0] == '>':
        header = i[1:].strip()
        if seq:
          kmerator.search_kmers(sequence.Sequence(header, seq))
          seq = ''
      else:
        seq += i.strip()
    kmerator.search_kmers(sequence.Sequence(header, seq))
    fh.close()

    def __init__(self):
      pass
