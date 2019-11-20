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

import kmerdb

class Kmerator:

  db = kmerdb.KmerDb()
  kmer_count = 0

  @staticmethod
  def calculate_kuid(kmer_seq):
    return hashlib.sha1(kmer_seq.encode()).hexdigest()

  def __init__(self, kmer_len, min_sequences_per_kmer, max_kmers_any_sequence):
    self.kmer_len = int(kmer_len)
    self.min_sequences_per_kmer = min_sequences_per_kmer
    self.max_kmers_any_sequence = max_kmers_any_sequence
    self.preselected_kmers = set()
    self.skip_kmers = set()
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
    if kuid not in self.skip_kmers:
      kmer = Kmerator.db.get_kmer(kuid, kmer_seq, seqname)
      location = Kmerator.db.new_kmer_location(kmer_start, seqname, kmer)
      neighbor_loc = self.find_neighbor_loc(location, kmer)
      if neighbor_loc and neighbor_loc.superkmer is None:
        skmer_idx = Kmerator.db.add_new_superkmer(neighbor_loc, location)
        if neighbor_loc.sequence != location.sequence:
          sys.exit("Error: Wrong superkmer locations on {} at {}: {}\t{}".format(seqname, kmer_start, neighbor_loc.sequence, location.sequence))
        #print("New superkmer:{} : {}\n\tsequence\t{}\n\tsloc\t{}\n\tstart\t{}\n\teloc\t{}\n\tstop\t{}".format(skmer_idx,kmer_seq,
                                                                                                          #neighbor_loc.sequence,
                                                                                                          #neighbor_loc.idx,
                                                                                                          #Kmerator.db.get_superkmer_start(skmer_idx).start,
                                                                                                          #location.idx,
                                                                                                          #Kmerator.db.get_superkmer_end(skmer_idx).end))
      if neighbor_loc and neighbor_loc.superkmer is not None:  # extend exisitng superkmer:
        skmer_idx = Kmerator.db.extend_superkmer(neighbor_loc.superkmer, location)
        #print("\tExtend superkmer {}\n\tfrom {}\n\tto {}\n\tby {}".format(skmer_idx, Kmerator.db.get_superkmer_start(skmer_idx).start,
                                                                           #Kmerator.db.get_superkmer_end(skmer_idx).end,
                                                                           #location.idx))
      kmer.add_location(location)
      # do the SWIGG min_alt_seqs and repeat_threshold_across check
      #self.preselect_kmer(kmer)
      Kmerator.kmer_count += 1

  def find_neighbor_loc(self, location, kmer):
    """
    Find the location of a kmer neighbor, i.e. the kmer preceding the current
    kmer by lmer_len.
    """
    if len(kmer.locations[location.sequence]) < 2:
      return None
    prev_kmer_location = Kmerator.db.get_kmer_sequence_location(location.sequence, (location.start - kmer.length() + 1))
    if not prev_kmer_location:
      return None
    #print(prev_kmer_location.start, prev_kmer_location.kuid, location.start, location.kuid)
    if prev_kmer_location.kuid == location.kuid:
      return prev_kmer_location
    return None

  def search_kmers(self, sequence):
    for i in range(sequence.length()-self.kmer_len):
      self.add_kmer(sequence.sequence[i:(i+self.kmer_len)], i, sequence.header)

  def show_kmers(self):
    print("Analyzed {} kmers. Kept {} distinct kmers in {} locations in {} sequences".format(Kmerator.kmer_count,
                                                                                             Kmerator.db.kmer_count(),
                                                                                             Kmerator.db.location_count(),
                                                                                             Kmerator.db.sequence_count()))
    Kmerator.db.show_kmers()

  def show_locations(self):
    print("Found {} kmers in {} locations".format(len(Kmerator.kmerdb), len(Kmerator.locdb)))
    for i in Kmerator.kmerdb:
      print(Kmerator.kmerdb[i].kuid, Kmerator.kmerdb[i].sequence, sep='\t')
      for j in Kmerator.kmerdb[i].locations:
        for k in Kmerator.kmerdb[i].locations[j]:
          print("\t\t\t", Kmerator.locdb[j][k].idx, Kmerator.locdb[j][k].start, Kmerator.locdb[j][k].end(), Kmerator.locdb[j][k].superkmer,  Kmerator.locdb[j][k].sequence, Kmerator.locdb[j][k], sep='\t')

  def show_superkmers(self):
    print("Found {} superkmers".format(Kmerator.db.superkmer_count()))
    Kmerator.db.show_superkmers()

  def preselect_kmer(self, kmer):
    if kmer.location_max_count > self.max_kmers_any_sequence:
      # This criteria won't change once it's reached and does not need any additional checks.
      self.skip_kmers.add(kmer.kuid)
      # rm kmer from preselected
      if kmer.kuid in self.preselected_kmers:
        self.preselected_kmers.remove(kmer.kuid)
      # rm kmer from database
      #print("Removing kmer from db: {}".format(kmer.kuid))
      Kmerator.db.remove_kmer(kmer)
    else:
      # Kmer number is less than requested for any sequence but maybe not found in the requested number of sequences.
      if len(kmer.locations) >= self.min_sequences_per_kmer:
        self.preselected_kmers.add(kmer.kuid)

  def filter(self, max_kmers_per_sequence):
    for i in self.preselected_kmers:
      if not self.hasMaxKmerPerSequence(Kmerator.kmerdb[i].locations, max_kmers_per_sequence):
        self.selected_kmers[i] = 0
      else:
        self.skip_kmers.add(i)
    locations_to_remove = []
    for i in self.skip_kmers:
      for j in Kmerator.kmerdb[i].locations:
        locations_to_remove += Kmerator.kmerdb[i].locations[j]
    print(len(locations_to_remove))
    print(locations_to_remove)

  def hasMaxKmerPerSequence(self, kmer_locations, max_kmers_per_sequence):
    for i in kmer_locations:
      if len(kmer_locations[i]) > max_kmers_per_sequence:
        return True
    return False
