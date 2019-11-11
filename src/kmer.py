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

    def __init__(self, start, kuid, sequence):
      self.start = start
      self.kuid = kuid
      self.sequence = sequence
      self.idx = Kmerator.KmerLocation.count
      Kmerator.KmerLocation.count += 1

  class Kmer:

    kmers = 0

    def __init__(self, kmer_seq, kuid):
      self.sequence = kmer_seq
      self.kuid = kuid
      self.count = 0
      self.locations = {}
      Kmerator.Kmer.kmers += 1

  def __init__(self, kmer_len, min_sequences_per_kmer, max_kmers_any_sequence):
    self.kmer_len = int(kmer_len)
    self.min_sequences_per_kmer = min_sequences_per_kmer
    self.max_kmers_any_sequence = max_kmers_any_sequence
    self.preselected_kmers = {}
    self.selected_kmers = {}

  def add_kmer(self, kmer_seq, kmer_start, seqname):
    kuid = Kmerator.calculate_kuid(kmer_seq)
    kmerlocation = Kmerator.KmerLocation(kmer_start, kuid, seqname)
    Kmerator.kmerlocdb[kmerlocation.idx] = kmerlocation
    if kuid  not in Kmerator.kmerdb:
      Kmerator.kmerdb[kuid] = Kmerator.Kmer(kmer_seq, kuid)
    Kmerator.kmerdb[kuid].count += 1
    if kmerlocation.sequence not in Kmerator.kmerdb[kuid].locations:
      Kmerator.kmerdb[kuid].locations[kmerlocation.sequence] = []
    Kmerator.kmerdb[kuid].locations[kmerlocation.sequence].append(kmerlocation.idx)
    # do the SWIGG min_alt_seqs and repeat_threshold_across check
    self.isInRepeatThresholds(kuid)

  def search_kmers(self, sequence):
    for i in range(sequence.length()-self.kmer_len):
      self.add_kmer(sequence.sequence[i:(i+self.kmer_len)], i, sequence.header)

  def show(self):
    print("Found {} kmers in {} locations".format(Kmerator.kmer_count(), Kmerator.location_count()))
    print("KmerUid", "kmer_seq", "locations", sep='\t')
    for i in Kmerator.kmerdb:
      print(Kmerator.kmerdb[i].kuid, Kmerator.kmerdb[i].sequence, sep='\t')
      for j in Kmerator.kmerdb[i].locations:
        for k in Kmerator.kmerdb[i].locations[j]:
          print("\t\t", Kmerator.kmerlocdb[k].sequence, Kmerator.kmerlocdb[k].start, Kmerator.kmerlocdb[k].idx, Kmerator.kmerlocdb[k], sep='\t')

  def isInRepeatThresholds(self, kuid):
    if (len(Kmerator.kmerdb[kuid].locations) >= self.min_sequences_per_kmer) and (len(Kmerator.kmerdb[kuid].locations) <= self.max_kmers_any_sequence):
      self.preselected_kmers[kuid] = 0
    else:
      if kuid in self.preselected_kmers:
        self.preselected_kmers.pop(kuid)

  def filter(self, max_kmers_per_sequence):
    for i in self.preselected_kmers:
      if self.hasMaxKmerPerSequence(Kmerator.kmerdb[i].locations, max_kmers_per_sequence):
        self.selected_kmers[i] = 0

  def hasMaxKmerPerSequence(self, kmer_locations, max_kmers_per_sequence):
    for i in kmer_locations:
      if len(kmer_locations[i]) > max_kmers_per_sequence:
        return False
    return True
