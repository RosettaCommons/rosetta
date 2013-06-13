#!/usr/bin/env python

def load_files():
	import json
	new=None
	try:
		new = json.load(  open('new/runtimes.yaml') )
	except:
		print 'Missing new/runtimes.yaml'
	ref=None
	try:
		ref = json.load( open('ref/runtimes.yaml') )
	except:
		print 'Missing new/runtimes.yaml'
	return new, ref


def compare_times():
	new, ref = load_files()	
	runtimes_compare = open('runtime_diffs.txt', 'w')
	format_string = '%39s %7.2f %7.2f %7.2f %7.2f\n'
	runtimes_compare.write('%39s %7s %7s %7s %7s\n' % ("TEST", "NEW", "REF", "DIFF", "D/R"))

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
			runtimes_compare.write( format_string % (key, n, r, diff, diff/r) )
			diffs.append(diff)
			new_values.append(n)
			ref_values.append(r)
	new_sum = sum(new_values)
	ref_sum = sum(ref_values)

	if len(new_values) > 0 and len(ref_values) > 0:
		new_avg = new_sum/float(len(new_values))
		ref_avg = ref_sum/float(len(ref_values))
	
		runtimes_compare.write( format_string % ('TOTAL', new_sum, ref_sum, new_sum - ref_sum, (new_sum-ref_sum)/ref_sum) )
		runtimes_compare.write( format_string % ('MEAN', new_avg, ref_avg, new_avg - ref_avg, (new_avg-ref_avg)/ref_avg ) )
	
def main():
	compare_times()

from sys import exit
if __name__ == "__main__":
	exit(main())
