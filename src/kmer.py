#!/usr/bin/env python3
#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


class Kmer:

  def __init__(self, kmer_seq, kuid):
    self.sequence = kmer_seq
    self.kuid = kuid
    self.count = 0
    self.location_max_count = 0
    self.locations = {}

  def length(self):
    return len(self.sequence)

  def add_location(self, location):
    self.locations[location.sequence].append(location.idx)
    self.location_max_count = max(self.location_max_count, len(self.locations[location.sequence]))
