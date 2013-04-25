#!/bin/env python


def load_files():
	import json
	from glob import glob
	new_results = glob('tests/*/output/.results.yaml')
	new={}
	for file in new_results:
		test = file.split('/')[1]
		try:
			new_dict = json.load(  open(file) ) 
			for key in new_dict:
				new[(test, key)]= new_dict[key]
		except:
			print 'Cannot load '+file
	ref_results = glob('tests/*/output/.old_results.yaml')
	ref={}
	for file in ref_results:
		test = file.split('/')[1]
		try:
			ref_dict = json.load(  open(file) ) 
			for key in ref_dict:
				ref[(test, key)]= ref_dict[key]
		except:
			print 'Cannot load '+ file
	return new, ref


def compare_times():
	new, ref = load_files()	
	diff_file = open('test_diffs.txt', 'w')
	format_string = '%39s %7.0f %7.0f %7.0f %7.2f %s\n'
	diff_file.write('%39s %7s %7s %7s %7s %s\n' % ("TEST", "NEW", "REF", "DIFF", "D/R", "TYPE"))

	from math import isnan, isinf
	new_max_mem=0
	ref_max_mem=0
	new_sum_sec=0
	ref_sum_sec=0
	test_count=0

	for key in sorted(new):# this key is a tuple
		n = new[key]
		if key in ref:
			test_count += 1
			r = ref[key]
			if key[1] == 'execution_time_s':
				new_sum_sec += n
				ref_sum_sec += r
			else:
				new_max_mem = max( n, new_max_mem )
				ref_max_mem = max( r, new_max_mem )
			if isnan(n+r) or isinf(n+r):
				continue # cannot compare
			diff = n - r
			d_diff = 0
			if r != 0: 
				d_diff = diff/float(r)
			diff_file.write( format_string % (key[0], n, r, diff, d_diff,key[1]) )
	new_avg_sec = new_sum_sec/float(test_count)
	ref_avg_sec = ref_sum_sec/float(test_count)
	sum_sec_diff = new_sum_sec - ref_sum_sec
	avg_sec_diff = new_avg_sec - ref_avg_sec
	max_mem_diff = new_max_mem - ref_max_mem
	d_sum_sec = 0
	if ref_sum_sec != 0: d_sum_sec = sum_sec_diff/float(ref_sum_sec)
	d_avg_sec = 0
	if ref_avg_sec != 0: d_avg_sec = avg_sec_diff/float(ref_avg_sec)
	d_max_mem = 0
	if ref_max_mem != 0: d_max_mem = max_mem_diff/float(ref_max_mem)
	
	diff_file.write( format_string % ('SUMMARY', new_sum_sec, ref_sum_sec, sum_sec_diff, d_sum_sec, 'Total Execution Time') )
	diff_file.write( format_string % ('SUMMARY', new_avg_sec, ref_avg_sec, avg_sec_diff, d_avg_sec, 'Average Execution Time' ) )
	diff_file.write( format_string % ('SUMMARY', new_max_mem, ref_max_mem, max_mem_diff, d_max_mem, 'Overall Max Memory' ) )
	
def main():
	compare_times()

from sys import exit
if __name__ == "__main__":
	exit(main())
