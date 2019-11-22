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
  """Manage kmer and their locations"""

  kmers = {}
  locations = {}
  sequence_locations = {}
  superkmers = {}
  superkmer_idx = 0

  @staticmethod
  def kmer_count():
    return len(KmerDb.kmers)

  @staticmethod
  def location_count():
    return len(KmerDb.locations)

  @staticmethod
  def sequence_count():
    return len(KmerDb.sequence_locations)

  @staticmethod
  def superkmer_count():
    return len(KmerDb.superkmers)

  def __init__(self):
    pass

  def get_kmer(self, kuid, kmer_seq, location):
    """
    Return a known kmer or create a new one. Kmers are indexed in KmerDb.kmers
    by their kuid (hash).

    :param hash kuid: uniqe hash for kmer sequence
    :param str kmer_seq: kmer sequence
    :param str location: sequence name on which the kmer was found
    :return: kmer
    :rtype: kmer.Kmer
    """
    if kuid not in KmerDb.kmers:
      KmerDb.kmers[kuid] = kmer.Kmer(kmer_seq, kuid)
    k = KmerDb.kmers[kuid]
    k.count += 1
    if location not in k.locations:
      k.locations[location] = []
    return k

  def new_location(self, kmer_start, kuid, seqname):
    """
    Return a new location for a kmer.

    :param int kmer_start: location start
    :param hash kuid: Kmer hash
    :param seqname str: sequence name on which the kmer was found
    :return: kmer location
    :rtype: kmerloc.Location
    """
    return kmerloc.KmerLocation(kmer_start, kuid, seqname)

  def add_kmer_location(self, location):
    """
    Add a new location for a kmer. All kmer locations are indexed in
    KmerDb.locations. Locations per sequence are referenced as
    start_coordinate:index in KmerDb.sequence_locations.
    KmerDb.locations indices are generated  automatically and not related to
    the sequence where they werer found on.

    :param int kmer_start: location start
    :param hash kuid: Kmer hash
    :param seqname str: sequence name on which the kmer was found
    :return: kmer location
    :rtype: kmerloc.Location
    """
    KmerDb.locations[location.idx] = location
    if location.sequence not in KmerDb.sequence_locations:
      KmerDb.sequence_locations[location.sequence] = {}
    KmerDb.sequence_locations[location.sequence][location.start] = location.idx


  def remove_locations(self, locations):
    """
    Remove locations. WIP

    :param list locations: location indices
    :return: kmer location
    :rtype: kmerloc.Location
    """
    print(locations)
    #return location

  def get_location_by_idx(self, idx):
    """
    Return location by its index from the location database.

    :param int idx: location idx
    :return: location, if found
    :rtype: kmerloc.Location or None
    """
    return KmerDb.locations.get(idx, None)

  def get_location_by_indices(self, indices):
    """
    Return location by its index from the location database.

    :param list indices: list with location indices
    :return: list of locations, if found
    :rtype: list of kmerloc.Location or None
    """
    locations = []
    for i in indices:
      if i in KmerDb.locations:
        locations.append(KmerDb.locations[i])
    return locations

  def get_kmer_sequence_locations(self, sequence_name, coords=None):
    """
    Get kmer locations for a given sequence. If no coordinate are given, all
    locations are returned.

    :param str sequence_name: name of sequence to look for locations
    :param list coords: location start coodinates
    :return: found locations
    :rtype: list of indices or None
    """
    if sequence_name in KmerDb.sequence_locations and coords is None:
      return [x for x in KmerDb.sequence_locations[sequence_name]]
    if sequence_name in KmerDb.sequence_locations:
      locations = []
      for i in coords:
        if i in KmerDb.sequence_locations[sequence_name]:
          locations.append(KmerDb.sequence_locations[sequence_name][i])
      return locations
    return None

  def get_kmer_sequence_location(self, sequence_name, coord):
    """
    Get a single kmer location for a given sequence.

    :param str sequence_name: name of sequence to look for locations
    :param int coord: location start coodinate
    :return: location if found or None
    :rtype: kmerloc.KmerLocation
    """
    if sequence_name in KmerDb.sequence_locations and (coord in KmerDb.sequence_locations[sequence_name]):
      return KmerDb.locations[KmerDb.sequence_locations[sequence_name][coord]]
    return None

  def count_kmers_on_sequence(self, sequence_name):
    """
    Return number of all found kmers on a given sequence

    :param str sequence_name: sequence name
    :return: number of kmer
    :rtype: int
    """
    if sequence_name in Kmerator.sequence_locations:
      return len(Kmerator.sequence_locations[sequence_name])
    return 0

  def add_new_superkmer(self, startloc, endloc):
    """
    Add new superkmer to the database KmerDb.superkmers. Superkmers are defined
    as multiple of the kmer len between two locations of the same kmer.
    """
    KmerDb.superkmer_idx += 1
    KmerDb.superkmers[KmerDb.superkmer_idx] = [startloc.idx, endloc.idx]
    startloc.superkmer = KmerDb.superkmer_idx
    endloc.superkmer = KmerDb.superkmer_idx
    return KmerDb.superkmer_idx

  def extend_superkmer(self, superkmer_idx, location):
    """Extend existing superkmer by adjusting the end location"""
    location.superkmer = superkmer_idx
    KmerDb.superkmers[superkmer_idx][1] = location.idx
    return superkmer_idx

  def get_superkmer_start(self, superkmer_idx):
    """Return the first location object of a superkmer"""
    return KmerDb.locations[KmerDb.superkmers[superkmer_idx][0]]

  def get_superkmer_end(self, superkmer_idx):
    """Return the last location object  of a superkmer"""
    return KmerDb.locations[KmerDb.superkmers[superkmer_idx][1]]

  def remove_superkmer(self, superkmer_idx):
    KmerDb.superkmers.pop(superkmer_idx, None)

  def remove_kmer(self, kmer):
    """Remove a kmer form the database. It has to be removed from all dicts.
    Rather inefficient so far."""
    kmer = KmerDb.kmers.pop(kmer.kuid)
    for i in kmer.locations:
      seqlocs = kmer.locations[i]
      #print("Removing {} / {} from sequence {}".format(kmer.kuid, kmer.sequence, i))
      for j in seqlocs:
        loc = KmerDb.locations.pop(j)
        #print("Removed loc {} for kmer {}: {} {} from locations".format(loc.idx, loc.kuid, loc.start, loc.end))
        KmerDb.sequence_locations[i].pop(loc.start)
        #print("Removed loc {} for kmer {} from sequence locations".format(loc.idx, loc.kuid, loc.start))
        if loc.superkmer is not None:
          #print("Removing superkmer {}".format(loc.superkmer))
          self.remove_superkmer(loc.superkmer)

  def show_kmers(self, kmers=None):
    "Print kmer to STDOUT."
    if kmers is not None:
      for i in KmerDb.kmers:
        if i in kmers:
          print(KmerDb.kmers[i].kuid, KmerDb.kmers[i].sequence, sep='\t')
          for j in KmerDb.kmers[i].locations:
            locs = KmerDb.kmers[i].locations[j]
            for k in self.get_location_by_indices(locs):
              print("\t", k.idx, k.start, k.end, k.superkmer, k.sequence,sep='\t')
    else:
      for i in KmerDb.kmers:
        print(KmerDb.kmers[i].kuid, KmerDb.kmers[i].sequence, sep='\t')
        for j in KmerDb.kmers[i].locations:
          locs = KmerDb.kmers[i].locations[j]
          for k in self.get_location_by_indices(locs):
            print("\t", k.idx, k.start, k.end, k.superkmer, k.sequence,sep='\t')

  def show_superkmers(self):
    for i in KmerDb.superkmers:
      skmer_start = self.get_superkmer_start(i)
      skmer_end = self.get_superkmer_end(i)
      kmer = KmerDb.kmers[skmer_start.kuid]
      print(i, kmer.kuid,  kmer.sequence, skmer_start.idx, skmer_start.start, skmer_end.idx, skmer_end.end)

  def dump_kmers(self, fout='kmer.dump', kmers=None):
    fh = open(fout, 'w')
    for i in KmerDb.kmers:
      if kmers and i in kmers:
        fh.write("{}\t{}\t{}".format(KmerDb.kmers[i].kuid, KmerDb.kmers[i].sequence, KmerDb.kmers[i].count))
        for j in KmerDb.kmers[i].locations:
          locs = KmerDb.kmers[i].locations[j]
          for k in self.get_location_by_indices(locs):
            fh.write("\t{}\t{}\t{}\t{}\t{}\n".format(k.idx, k.start, k.end, k.superkmer, k.sequence))
      else:
        fh.write("{}\t{}\t{}\n".format(KmerDb.kmers[i].kuid, KmerDb.kmers[i].sequence, KmerDb.kmers[i].count))
        for j in KmerDb.kmers[i].locations:
          locs = KmerDb.kmers[i].locations[j]
          for k in self.get_location_by_indices(locs):
            fh.write("\t{}\t{}\t{}\t{}\t{}\n".format(k.idx, k.start, k.end, k.superkmer, k.sequence))
    fh.close()
