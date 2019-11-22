#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


class KmerLocation:

  idx = 0

  def __init__(self, sequence, start, kmer):
    self.start = start
    self.end = start + kmer.length() -1
    self.kuid = kmer.kuid
    self.sequence = sequence
    self.idx = KmerLocation.idx
    self.superkmer = None
    KmerLocation.idx += 1
