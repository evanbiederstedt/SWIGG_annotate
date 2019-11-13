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

  @staticmethod
  def calculate_kuid(kmer_seq):
    return hashlib.sha1(kmer_seq.encode()).hexdigest()

  @staticmethod
  def kmer_count():
    return Kmerator.Kmer.kmers

  @staticmethod
  def location_count():
    return Kmerator.KmerLocation.count

  class KmerLocation:

    count = 0
    idx = 0

    def __init__(self, start, kuid, sequence):
      self.start = start
      self.kuid = kuid
      self.sequence = sequence
      self.idx = Kmerator.KmerLocation.idx
      Kmerator.KmerLocation.idx += 1
      Kmerator.KmerLocation.count += 1

    def end(self):
      return self.start + Kmerator.kmerdb[self.kuid].length() - 1

  class Kmer:

    kmers = 0

    def __init__(self, kmer_seq, kuid):
      self.sequence = kmer_seq
      self.kuid = kuid
      self.count = 0
      self.location_max = 0
      self.locations = {}
      self.seqrep = 1
      Kmerator.Kmer.kmers += 1

    def length(self):
      return self.seqrep * len(self.sequence)

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
      kmer = Kmerator.kmerdb.get(kuid, None)
      if not kmer:
        kmer = Kmerator.Kmer(kmer_seq, kuid)
      kmerlocation = Kmerator.KmerLocation(kmer_start, kuid, seqname)
      lh_neighbor_loc = self.find_neighbor_location(kmer, kmerlocation)
      # if kmer has neighbor, merge and proceed with merged kmer
      if lh_neighbor_loc:
        kmer = self.merge_neighbor_kmer(kmer, kmerlocation, lh_neighbor_loc)
        kmerlocation = Kmerator.KmerLocation(lh_neighbor_loc.start, kmer.kuid, lh_neighbor_loc.sequence)
        Kmerator.KmerLocation.count -= 1
      if kmerlocation.sequence not in kmer.locations:
        kmer.locations[kmerlocation.sequence] = []
      kmer.locations[kmerlocation.sequence].append(kmerlocation.idx)
      kmer.count += 1
      if kmer.kuid not in Kmerator.kmerdb:
        Kmerator.kmerdb[kmer.kuid] = kmer
      Kmerator.kmerlocdb[kmerlocation.idx] = kmerlocation
      kmer.location_max = max(kmer.location_max, len(kmer.locations[kmerlocation.sequence]))
      # do the SWIGG min_alt_seqs and repeat_threshold_across check
      self.preselect_kmer(kmer)

  def find_neighbor_location(self, kmer, location):
    # If we haven't seen the kmer on
    if not location.sequence in kmer.locations:
      return None
    if kmer.locations[location.sequence]:
      print(kmer.locations[location.sequence])
      # Fetch latest added location of kmer on sequence
      prev_kmer_loc =  Kmerator.kmerlocdb[kmer.locations[location.sequence][-1]]
      print("N-Test:\nkmer0: {}\t{}\t{}\t{}\nkmer1: {}\t{}\t{}\t{}".format(prev_kmer_loc.kuid,
                                                                          prev_kmer_loc.start,
                                                                          prev_kmer_loc.end(),
                                                                          prev_kmer_loc.sequence,
                                                                          location.kuid,
                                                                          location.start,
                                                                          location.end(),
                                                                          location.sequence))
      if location.start - prev_kmer_loc.end() == 1:
        return prev_kmer_loc
    return None

  def merge_neighbor_kmer(self, kmer, location, lh_neighbor_loc):
    merged_kmer = Kmerator.Kmer(kmer.sequence, Kmerator.calculate_kuid(''.join([kmer.sequence, kmer.sequence])))
    merged_kmer.seqrep += 1
    if Kmerator.kmerdb[kmer.kuid].count - 1 == 0:
      Kmerator.kmerdb.pop(kmer.kuid)
      print("Rm kmer:",  kmer.kuid)
      self.preselected_kmers.pop(kmer.kuid)
    else:
      Kmerator.kmerdb[kmer.kuid].count -= 1
    Kmerator.Kmer.kmers -= 1
    print("Merged: {} into {}".format(kmer.kuid, merged_kmer.kuid))
    return merged_kmer

  def search_kmers(self, sequence):
    for i in range(sequence.length()-self.kmer_len):
      self.add_kmer(sequence.sequence[i:(i+self.kmer_len)], i, sequence.header)

  def show(self):
    print("Found {} kmers in {} locations".format(Kmerator.kmer_count(), Kmerator.location_count()))
    for i in Kmerator.kmerdb:
      print(Kmerator.kmerdb[i].kuid, Kmerator.kmerdb[i].sequence, sep='\t')
      for j in Kmerator.kmerdb[i].locations:
        for k in Kmerator.kmerdb[i].locations[j]:
          print("\t\t\t", Kmerator.kmerlocdb[k].sequence, Kmerator.kmerlocdb[k].start, Kmerator.kmerlocdb[k].idx, Kmerator.kmerlocdb[k], sep='\t')

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
