#  -------------------------------------------------------------------------------
#  \author Jan P Buchmann <jan.buchmann@sydney.edu.au>
#  \copyright 2019 The University of Sydney
#  \description
#  -------------------------------------------------------------------------------


import io
import os
import sys

class Sequence:

  def __init__(self, header, sequence, metadata=None):
    self.header = header
    self.sequence = sequence
    self.metadata = metadata

  def length(self):
    return len(self.sequence)
