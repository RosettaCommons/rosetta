#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

# Utility script to parse and extract data from score files in JSON format
# @author Luki Goldschmidt <lugo@uw.edu>

import sys
import os
import json
from collections import OrderedDict
from optparse import OptionParser

class ScoreFile:
  def __init__(self, filename):
    self.filename = filename
    self.decoys = []
    self.decoy_field_name = "decoy"

    if filename == "-":
      lines = file(sys.stdin).readlines()
    else:
      lines = file(filename).readlines()

    for line in lines:
      try:
        o = json.loads(line)
        self.decoys.append(o)
      except ValueError:
        print >> sys.stderr, "Failed to parse JSON object; skipping line:\n", line

  def getDecoyCount(self):
    return len(self.decoys)

  def getDecoyNames(self):
    return [ str(r[self.decoy_field_name]) for r in self.decoys ]

  def getScoreTermNames(self):
    r = []
    for rec in self.decoys:
      for k in rec.keys():
        if k != self.decoy_field_name and not k in r:
          r.append(str(k))
    r.sort()
    return r

  def getScoreTerms(self, scoreterms=""):
    if type(scoreterms) == str:
      if scoreterms in ["", "*"]:
        scoreterms = self.getScoreTermNames()
        scoreterms.sort()
      else:
        scoreterms = scoreterms.split(",")
    r = OrderedDict()
    for rec in self.decoys:
      scores = {}
      for scoreterm in scoreterms:
        scores[scoreterm] = rec.get(scoreterm)
      r[ str(rec[self.decoy_field_name]) ] = scores
    return r

  def getScoreTerm(self, scoreterm):
    r = {}
    for rec in self.decoys:
      r[ str(rec[self.decoy_field_name]) ] = rec.get(scoreterm)
    return r

  def getStats(self, scoreterms=""):
    scores = self.getScoreTerms(scoreterms)
    stats = None

    # Collect data
    for decoy_name in scores:
      decoy = scores[decoy_name]
      if stats == None:
        columns = decoy.keys()
        stats = { k: [] for k in columns }

      for term in decoy:
        score = decoy[term]
        if not score is None and type(score) == float:
          stats[term].append(score)

    # Compute Stats
    calc_stats = OrderedDict()

    def median(s):
      s.sort()
      return s[ int(len(s)/2) ]

    def stddev(s):
      mean = sum(s) / len(s)
      r = 0.0
      for x in s:
        r += pow(x - mean, 2)
      return pow(r / len(s), 0.5)

    for column in columns:
      s = stats[column]
      calc_stats[column] = {
        'n':      len(s),
        'mean':   sum(s) / len(s) if len(s) > 0 else None,
        'median': median(s)       if len(s) > 0 else None,
        'stddev': stddev(s)       if len(s) > 0 else None,
        'min':    min(s)          if len(s) > 0 else None,
        'max':    max(s)          if len(s) > 0 else None
      }

    return calc_stats

########################################################################

def output_legacy(out, prefix):
  # Format floats
  for i in range(0, len(out)):
    for j in range(0, len(out[i])):
      t = out[i][j]
      value = t[1]
      if type(value) == float:
        value = "%.03f" % value
        out[i][j] = (t[0], value)

  # Compute column widths
  widths = {}
  for decoy_terms in out:
    for decoy_term in decoy_terms:
      name, value = decoy_term
      if not widths.has_key(name):
        widths[name] = len(name)
      l = len(str(value))
      if widths[name] < l:
        widths[name] = l

  if widths.has_key('decoy'):
    widths['decoy'] = -widths['decoy']

  # Header
  sys.stdout.write(prefix)
  for decoy_term in out[0]:
    name, value = decoy_term
    sys.stdout.write("%*s " % (widths[name], name))
  sys.stdout.write("\n")

  # Scores
  for decoy_terms in out:
    sys.stdout.write(prefix)
    for decoy_term in decoy_terms:
      name, value = decoy_term
      sys.stdout.write("%*s " % (widths[name], value))
    sys.stdout.write("\n")

def output_tab(out, prefix):
  # Header
  names = [ decoy_term[0] for decoy_term in out[0] ]
  sys.stdout.write(prefix + "\t".join(names) + "\n")

  # Scores
  for decoy_terms in out:
    values = [ str(decoy_term[1]) for decoy_term in decoy_terms ]
    sys.stdout.write(prefix + "\t".join(values) + "\n")

def output_CSV(out, prefix):
  # Header
  names = [ decoy_term[0] for decoy_term in out[0] ]
  sys.stdout.write(prefix + ",".join(names) + "\n")

  # Scores
  for decoy_terms in out:
    values = []
    for decoy_term in decoy_terms:
      value = decoy_term[1]
      if type(value) in [str, unicode] and "," in value:
        value = '"%s"' % value
      values.append(str(value))
    sys.stdout.write(prefix + ",".join(values) + "\n")

def printVerbose(s):
  if options.verbose:
    print >> sys.stderr, s

########################################################################

def main(argv):
  parser = OptionParser(usage="usage: %prog [OPTIONS] scorefile [scorefile...]\n" +
    "This utility parses and extracts data from score files in JSON format")
  parser.set_description(main.__doc__)

  parser.add_option("-n", "--decoynames",
    action="store_true",
    default=False,
    help="List decoy names",
  )

  parser.add_option("-t", "--scoreterms",
    action="store_true",
    default=False,
    help="List score term names",
  )

  parser.add_option("-S", "--summary",
    action="store_true",
    default=False,
    help="Compute stats summarizing data",
  )

  parser.add_option("-s", "--scores",
    default="",
    help="Comma separated list of score terms to extract",
  )

  parser.add_option("-o", "--output",
    default="legacy",
    help="Output format: legacy, tab, CSV",
  )

  parser.add_option("-p", "--prefix",
    default=None,
    help="Line prefix in score output (legacy format)",
  )

  parser.add_option("-v", "--verbose",
    action="store_true",
    default=False,
    help="Verbose output",
  )

  global options
  (options, remaining_args) = parser.parse_args(args=argv)

  if len(remaining_args) < 1:
    parser.parse_args(["--help"])
    sys.exit(1)

  for filename in remaining_args:
    if filename != "":
      printVerbose("    Scorefile: %s" % filename)

    if filename != "-" and not os.path.isfile(filename):
      print >> sys.stderr, "File not found:", filename
      continue

    sf = ScoreFile(filename)
    if options.verbose:
      printVerbose("       Decoys: %d" % sf.getDecoyCount())
      printVerbose("  Score terms: %s" % ", ".join(sf.getScoreTermNames()))
      printVerbose("")

    ### Info handlers
    if options.decoynames:
      print "\n".join( sf.getDecoyNames() )
      continue

    if options.scoreterms:
      print "\n".join( sf.getScoreTermNames() )
      continue

    ### Stats summary
    if options.summary:
      stats = sf.getStats(options.scores)
      max_width = max( [ len(x) for x in stats.keys() ] )
      fmt = "%*s:  %4s  %10s  %10s  %10s  %10s  %10s"
      print fmt % (max_width, "TERM", "n", "Min", "Max", "Mean", "Median", "StdDev")

      for column in stats:
        v = stats[column]
        for f in ['min','max','mean','median','stddev']:
          if not v[f] == None:
            v[f] = "%.3f" % v[f]
        print fmt % (max_width, column, v['n'], v['min'], v['max'], v['mean'], v['median'], v['stddev'])
      continue

    ### Default score list handler
    scores = sf.getScoreTerms(options.scores)
    out = []
    for decoy_name in scores:
      terms = scores[decoy_name]
      decoy_scores = [("decoy", decoy_name)]
      for term in terms:
        decoy_scores.append((term, terms[term]))
      out.append(decoy_scores)

    ### Output
    if options.output == "legacy":
      output_legacy(out, "SCORE: " if options.prefix is None else options.prefix + ": ")
    elif options.output == "tab":
      output_tab(out, "" if options.prefix is None else options.prefix + "\t")
    elif options.output == "CSV" or options.output == "csv":
      output_CSV(out, "" if options.prefix is None else options.prefix)

########################################################################

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
