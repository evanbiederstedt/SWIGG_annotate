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
import networkx

import kmerdb
from graph import node

class Kmerator:
  """Searching, filtering, and inspecting kmers in sequences."""

  db = kmerdb.KmerDb()
  kmer_count = 0

  @staticmethod
  def calculate_kuid(kmer_seq):
    """Calculate has for kmer sequence"""
    return hashlib.sha1(kmer_seq.encode()).hexdigest()

  def __init__(self, kmer_len, min_sequences_per_kmer, max_kmers_any_sequence):
    self.kmer_len = int(kmer_len)
    self.min_sequences_per_kmer = min_sequences_per_kmer
    self.max_kmers_any_sequence = max_kmers_any_sequence
    self.preselected_kmers = set()
    self.skip_kmers = set()
    self.selected_kmers = {}
    self.orderd_locations = {}

  def add_kmer(self, kmer_seq, kmer_start, seqname):
    """
    Check kmer and its locations for neighbors and do prselection based on kmer
    selection criteria. It tests if  kmer already exists, prepares its loation,
    tests for a possible kmer neighbor. If a direct kmer neighbor is found,
    creates or updates the corresponding superkmer. If the kmer passes
    prefiltering, it's added to the database.

    :param str kmer_seq: kmer sequence
    :param int kmer_start: kmer start coordinate
    :param str seqname: sequence name on which the kmer was found
    :return: kmer location or None
    :rtype: kmerloc.KmerLocation
    """
    kuid = Kmerator.calculate_kuid(kmer_seq)
    if kuid in self.skip_kmers:
      Kmerator.db.kmers[kuid].count += 1
      return None
    kmer = Kmerator.db.get_kmer(kuid, kmer_seq, seqname)
    location = Kmerator.db.new_location(seqname, kmer_start, kmer)
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
    Kmerator.kmer_count += 1
    # do the SWIGG min_alt_seqs and repeat_threshold_across check
    if self.preselect_kmer(kmer):
      kmer.add_location(location)
      Kmerator.db.add_kmer_location(location)
      return location
    return None

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
    if prev_kmer_location.kuid == location.kuid:
      return prev_kmer_location
    return None

  def search_kmers(self, sequence):
    """Scan sequence for kmers and store found kmers in order of occurence"""
    if sequence.header not in self.orderd_locations:
      self.orderd_locations[sequence.header] = []
    for i in range(sequence.length()-self.kmer_len):
      location = self.add_kmer(sequence.sequence[i:(i+self.kmer_len)], i, sequence.header)
      if location:
        self.orderd_locations[location.sequence].append(location)

  def find_anchor_kmers(self, fout='kmers.xml'):
    """
    Scan found kmers and assemble nodes for graph of anchor nodes. Add metadata
    for each node from the corresspoding kmer
    """
    nodes = {}
    G = networkx.DiGraph()
    for i in self.orderd_locations:
      loc = self.orderd_locations[i]
      prev = loc[0]
      parent = None
      for j in loc[1:]:
        if j.kuid in self.preselected_kmers and prev.kuid in self.preselected_kmers:
          if j.start - prev.start > self.kmer_len:
            if prev.kuid not in nodes:
              nodes[prev.kuid] = node.Node(prev.kuid, Kmerator.db.kmers[prev.kuid].sequence, Kmerator.db.kmers[prev.kuid].count, Kmerator.db.kmers[prev.kuid].sequence_count())
            n = nodes[prev.kuid]
            G.add_node(n.seq, count=n.count, sequences=n.sequence_count)
            if parent is not None:
              if parent.kuid != n.kuid:
                n.incoming.add(parent.kuid)
                G.add_edge(n.seq, parent.seq)
                if parent.kuid in nodes:
                  nodes[parent.kuid].out.add(n.kuid)
            parent = n
            #print(prev.idx, prev.kuid, prev.start, prev.end, Kmerator.db.kmers[prev.kuid].sequence, "<- this kmer")
            #print(j.idx, j.kuid, j.start, j.end, Kmerator.db.kmers[j.kuid].sequence)
            #print("++++++++++++++++++++++++++")
          prev = j
      #for i in nodes:
        #nodes[i].dump()
        #print("~~~~~~~~")
      #print("--------------------------------")
      networkx.write_graphml(G, fout)

  def show_kmers(self, kmers=None):
    if kmers is not None:
      print("Dumping {} kmers from database.".format(len(kmers)))
    else:
      print("Dumping all {} kmers in database.".format(Kmerator.db.kmer_count()))
    Kmerator.db.show_kmers(kmers)

  def show_locations(self):
    print("Found {} kmers in {} locations".format(len(Kmerator.kmerdb), len(Kmerator.locdb)))
    for i in Kmerator.kmerdb:
      print(Kmerator.kmerdb[i].kuid, Kmerator.kmerdb[i].sequence, sep='\t')
      for j in Kmerator.kmerdb[i].locations:
        for k in Kmerator.kmerdb[i].locations[j]:
          print("\t\t\t", Kmerator.locdb[j][k].idx, Kmerator.locdb[j][k].start, Kmerator.locdb[j][k].end(), Kmerator.locdb[j][k].superkmer,  Kmerator.locdb[j][k].sequence, Kmerator.locdb[j][k], sep='\t')

  def show_superkmers(self):
    """Show identified superkmers"""
    print("Found {} superkmers".format(Kmerator.db.superkmer_count()))
    Kmerator.db.show_superkmers()

  def preselect_kmer(self, kmer):
    """
    Test kmer for given thresholds and add to ignore list when thresholds have
    been crossed.
    """
    if self.trigger_repeat_threshold_across(kmer):#if triggers, don't preselect
      return False
    # Kmer number is less than requested for any sequence but maybe not found in the requested number of sequences, yet.
    if len(kmer.locations) >= self.min_sequences_per_kmer:
      self.preselected_kmers.add(kmer.kuid)
    return True

  def trigger_repeat_threshold_across(self, kmer):
    """
    Test for kemr repeat across any sequence. Once the threshold is reached, the
    kmer can be ignored and does not need any additional checks. If the kmer
    has been already preselcted, remove it from the selection"""
    if kmer.location_max_count > self.max_kmers_any_sequence:
      self.skip_kmers.add(kmer.kuid)
      if kmer.kuid in self.preselected_kmers:
        self.preselected_kmers.remove(kmer.kuid)
      #print("Trigger -ra: {} / {} due to -ra {}. Currently: {}".format(kmer.kuid, kmer.sequence, self.max_kmers_any_sequence, kmer.location_max_count, kmer.count))
      return True
    return False

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

  def dump_kmers(self, fout='kmer.dump'):
    """Dump kmers for debugging"""
    Kmerator.db.dump_kmers(fout)
