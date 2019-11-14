#!/usr/bin/env python3
#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


import io
import os
import sys
import hashlib
import random
import cProfile

class Kmerator:

  kmerdb = {}
  kmerlocdb = {}
  seqdb = {}
  superkmer = {}
  locdb = {}

  @staticmethod
  def calculate_kuid(kmer_seq):
    return hashlib.sha1(kmer_seq.encode()).hexdigest()

  @staticmethod
  def kmer_count():
    return len(Kmerator.kmerdb)

  @staticmethod
  def location_count():
    return Kmerator.KmerLocation.count

  class KmerLocation:

    idx = 0

    def __init__(self, start, kuid, sequence):
      self.start = start
      self.kuid = kuid
      self.sequence = sequence
      self.idx = Kmerator.KmerLocation.idx
      self.superkmer = None
      Kmerator.KmerLocation.idx += 1

    def end(self):
      return self.start + Kmerator.kmerdb[self.kuid].length() - 1

  class Kmer:

    def __init__(self, kmer_seq, kuid):
      self.sequence = kmer_seq
      self.kuid = kuid
      self.count = 0
      self.location_max = 0
      self.locations = {}

    def length(self):
      return len(self.sequence)

  def __init__(self, kmer_len, min_sequences_per_kmer, max_kmers_any_sequence):
    self.kmer_len = int(kmer_len)
    self.min_sequences_per_kmer = min_sequences_per_kmer
    self.max_kmers_any_sequence = max_kmers_any_sequence
    self.preselected_kmers = {}
    self.skipmap = {}
    self.selected_kmers = {}

  def add_kmer(self, kmer_seq, kmer_start, seqname):
    """
    - Test if kmer already exists
    - Prepare location
    - Test if kmer has direct neighbor
       Merge if yes
    - Add kmer and location to db
    """
    kuid = Kmerator.calculate_kuid(kmer_seq)
    if kuid not in self.skipmap:
      kmer = self.get_kmer(kuid, kmer_seq)
      location = self.get_kmer_location(kmer_start, kmer, seqname)
      neighbor_loc = self.find_neighbor_loc(location, kmer)
      if neighbor_loc:
        print("Neighbor kmer", neighbor_loc.start, neighbor_loc.end(), neighbor_loc.kuid, location.start, location.end(), location.kuid)
        self.adjust_superkmer(kmer, location, neighbor_loc)
      Kmerator.kmerlocdb[location.idx] = location
      kmer.locations[location.sequence].append(location.idx)
      kmer.location_max = max(kmer.location_max, len(kmer.locations[location.sequence]))
      if location.superkmer:
        print(location.superkmer, location.end())
      # do the SWIGG min_alt_seqs and repeat_threshold_across check
      #self.preselect_kmer(kmer)

  def get_kmer(self, kuid, kmer_seq):
    if kuid not in Kmerator.kmerdb:
      Kmerator.kmerdb[kuid] = Kmerator.Kmer(kmer_seq, kuid)
    Kmerator.kmerdb[kuid].count += 1
    return Kmerator.kmerdb[kuid]

  def get_kmer_location(self, kmer_start, kmer, seqname):
    kmerlocation = Kmerator.KmerLocation(kmer_start, kmer.kuid, seqname)
    if kmerlocation.sequence not in Kmerator.locdb:
      Kmerator.locdb[kmerlocation.sequence] = {}
    Kmerator.locdb[kmerlocation.sequence][kmerlocation.start] = kmerlocation
    if kmerlocation.sequence not in kmer.locations:
      kmer.locations[kmerlocation.sequence] = []
    return kmerlocation

  def adjust_superkmer(self, kmer, location, neighbor_loc):
    """Add new neigbbor (start new superkmer) or expand existing """
    if neighbor_loc.superkmer is None:  # start superkmer with neighbor_loc as start
      neighbor_loc.superkmer = len(Kmerator.superkmer)
      location.superkmer = neighbor_loc.superkmer
      Kmerator.superkmer[location.superkmer] = [neighbor_loc.start, location.end()]
      print("New superkmer: from {} to {}".format(Kmerator.superkmer[location.superkmer][0], Kmerator.superkmer[location.superkmer][1]))
    else: #neighbor_loc is already part of a superkmer
      print("Extend superkmer {} with {} to {}".format(neighbor_loc.superkmer, neighbor_loc.idx, location.end()))
      location.superkmer = neighbor_loc.superkmer
      Kmerator.superkmer[location.superkmer][1] = location.end()
      print("Coords for {}: {} to {}".format(neighbor_loc.superkmer, Kmerator.superkmer[neighbor_loc.superkmer][0], Kmerator.superkmer[location.superkmer][1]))

  def find_neighbor_loc(self, location, kmer):
    # test if kmer is part of a superkmer, i.e. more than two neighbors
    # find first neighbor
    if len(Kmerator.locdb[location.sequence]) < 2:
      return None
    if location.start % self.kmer_len != 0:
      return None
    if Kmerator.locdb[location.sequence][(location.start-self.kmer_len)].kuid == location.kuid:
      return Kmerator.locdb[location.sequence][(location.start-self.kmer_len)]
    return None

  def search_kmers(self, sequence):
    for i in range(sequence.length()-self.kmer_len):
      self.add_kmer(sequence.sequence[i:(i+self.kmer_len)], i, sequence.header)

  def show(self):
    print("Found {} kmers in {} locations".format(len(Kmerator.kmerdb), len(Kmerator.locdb)))
    for i in Kmerator.kmerdb:
      print(Kmerator.kmerdb[i].kuid, Kmerator.kmerdb[i].sequence, sep='\t')
      for j in Kmerator.kmerdb[i].locations:
        for k in Kmerator.kmerdb[i].locations[j]:
          print("\t\t\t", Kmerator.kmerlocdb[k].idx, Kmerator.kmerlocdb[k].start, Kmerator.kmerlocdb[k].end(), Kmerator.kmerlocdb[k].superkmer,  Kmerator.kmerlocdb[k].sequence, Kmerator.kmerlocdb[k], sep='\t')

  def preselect_kmer(self, kmer):
    if kmer.location_max > self.max_kmers_any_sequence:
      # This criteria won't change once it's reached and does not need any additional checks.
      self.skipmap[kmer.kuid] = 0
      # rm kmers from db
      if kmer.kuid in self.preselected_kmers:
        self.preselected_kmers.pop(kmer.kuid)
    else:
      # Kmer number is less than requested for any sequence but maybe not found in the requested number of sequences.
      if len(kmer.locations) >= self.min_sequences_per_kmer:
        self.preselected_kmers[kmer.kuid] = 0

  def filter(self, max_kmers_per_sequence):
    for i in self.preselected_kmers:
      if self.hasMaxKmerPerSequence(Kmerator.kmerdb[i].locations, max_kmers_per_sequence):
        self.selected_kmers[i] = 0

  def hasMaxKmerPerSequence(self, kmer_locations, max_kmers_per_sequence):
    for i in kmer_locations:
      if len(kmer_locations[i]) > max_kmers_per_sequence:
        return False
    return True
