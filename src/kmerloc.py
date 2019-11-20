#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


class KmerLocation:

  idx = 0

  def __init__(self, start, sequence, kmer):
    self.start = start
    self.end = kmer.length()
    self.kuid = kmer.kuid
    self.sequence = sequence
    self.idx = KmerLocation.idx
    self.superkmer = None
    KmerLocation.idx += 1

  #def end(self):
    #return self.start + Kmerator.kmerdb[self.kuid].length() - 1
