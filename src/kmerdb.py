#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


import io
import os
import sys

import kmer
import kmerloc

class KmerDb:

  kmers = {}
  locations = {}
  sequence_locations = {}

  def __init__(self):
    pass

  def get_kmer(self, kuid, kmer_seq, location):
    if kuid not in KmerDb.kmers:
      KmerDb.kmers[kuid] = kmer.Kmer(kmer_seq, kuid)
    k = KmerDb.kmers[kuid]
    k.count += 1
    if location.sequence not in k.locations:
      k.locations[location.sequence] = []
    return k


  def get_kmer_location(self, kmer_start, kuid, seqname):
    location = kmerloc.KmerLocation(kmer_start, seqname, KmerDb.kmers[kuid])
    Kmerator.locations[location.idx] = location
    if location.sequence not in KmerDb.sequence_locations:
      Kmerator.sequence_locations[location.sequence] = {}
    Kmerator.sequence_locations[location.sequence][kmerlocation.start] = location.idx
    return location

  def get_location_by_idx(self, idx):
    return KmerDb.locations.get(idx, None)

  def get_sequence_locations(self, sequence_name):
    if sequence_name in Kmerator.sequence_locations:
      return [x for x in Kmerator.sequence_locations[sequence_name]]
    return None
