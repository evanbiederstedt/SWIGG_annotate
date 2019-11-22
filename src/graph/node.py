#-------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#-------------------------------------------------------------------------------

class Node:

  def __init__(self, kuid, seq, count, sequence_count):
    self.kuid = kuid
    self.count = count
    self.seq = seq
    self.sequence_count = sequence_count
    self.incoming = set()
    self.out = set()

  def dump(self):
    return {'kuid': self.kuid,
            'seq': self.seq,
            'seq_count' : self.sequence_count,
            'in' : [x for x in self.incoming],
            'out' : [x for x in self.out]}
