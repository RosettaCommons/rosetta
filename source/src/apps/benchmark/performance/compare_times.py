#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

def load_files():
	import json
	new=None
	from os.path import exists
	try:
		new = json.load(  open('_performance_') )
	except:
		print 'Missing "_performance_" file'
	ref=None
	try:
		ref = json.load( open('_old_performance_') )
	except:
		print 'Missing "old_performance_" file'
	return new, ref


def compare_times():
	new, ref = load_files()
	runtimes_compare = open('runtime_diffs.txt', 'w')
	format_string = '%7.1f %7.1f %7.1f %6.2f %s\n'
	runtimes_compare.write('%7s %7s %7s %6s %s\n' % ("NEW", "REF", "DIFF", "D/R", "TEST"))

	from math import isnan, isinf
	diffs=[]
	new_values=[]
	ref_values=[]
	for key in set(new.keys() + ref.keys()):
		if key in new and key in ref:
			n = new[key]
			r = ref[key]
			if isnan(n+r) or isinf(n+r):
				continue # cannot compare
			diff = n -r
			runtimes_compare.write( format_string % (n, r, diff, diff/r, key) )
			diffs.append(diff)
			new_values.append(n)
			ref_values.append(r)
	new_sum = sum(new_values)
	ref_sum = sum(ref_values)
	new_avg = new_sum/float(len(new_values))
	ref_avg = ref_sum/float(len(ref_values))

	runtimes_compare.write( format_string % (new_sum, ref_sum, new_sum - ref_sum, (new_sum-ref_sum)/ref_sum, 'TOTAL' ))
	runtimes_compare.write( format_string % (new_avg, ref_avg, new_avg - ref_avg, (new_avg-ref_avg)/ref_avg, 'MEAN' ))

def main():
	compare_times()

from sys import exit
if __name__ == "__main__":
	exit(main())
