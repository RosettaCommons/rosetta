#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rosetta_cloud.py
## @brief  Provide API and command-line access to RosettaCloud
## @author Sergey Lyskov

import os, os.path, sys, zlib, json, copy, subprocess, shutil, urllib.parse

import time as time_module

from argparse import ArgumentParser
from configparser import ConfigParser

import requests
#import pdfkit
#import weasyprint
#from weasyprint import HTML
from datetime import datetime

################################################
def Sleep(time_ : int, message : str, dict_={}):
	''' Fancy sleep function '''
	len_ = 0
	for i in range(time_, 0, -1):
		#print "Waiting for a new revision:%s... Sleeping...%d	 \r" % (sc.revision, i),
		msg = message.format(dict_, time_left=i)
		print(msg, end='')
		len_ = max(len_, len(msg))
		sys.stdout.flush()
		time_module.sleep(5)

	print( ' '*len_ + '\r', end='') # erasing sleep message

################################################
def http_request(function, name, until_successes=True, **args):
	''' request URL using function method and return it content. Automatically initiate retry on failures.
	'''
	while True:
		try:
			r = function(_server_url_ + '/api/v1' + name, auth=(Config.get('main', 'user'), Config.get('main', 'password')), timeout=64, **args)
			if r.status_code == 200  or  (not until_successes): return r
			print('Request for {name} with extra arguments {args}\nReturned error {r.status_code} with message: {r.text}. Going to sleep then retry...'.format_map(locals()))
		except requests.exceptions.ConnectionError as e:
			print('Request for {name} with extra arguments {args}\nRaised exception {e}. Going to sleep then retry...'.format_map(locals()))

		Sleep(20, "Sleeping {time_left}... then I will retry...   \r" )

################################################
def get_url (name, **args):   return http_request(requests.get,	name, **args)
def post_url(name, **args):   return http_request(requests.post,   name, **args)
def delete_url(name, **args): return http_request(requests.delete, name, **args)

################################################
def get_test_ids():

	# download file
	os.system( "wget \'https://graylab.jhu.edu/Sergey/for.Julia/julia.txt\' --no-check-certificate -O tests.txt" )

	# read columns
	testid = subprocess.getoutput("awk '{print $1}' tests.txt").splitlines()
	name = subprocess.getoutput("awk '{print $2}' tests.txt").splitlines()
	teststatus = subprocess.getoutput("awk '{print $3}' tests.txt").splitlines()
	revision = subprocess.getoutput("awk '{print $4}' tests.txt").splitlines()
	testdir = subprocess.getoutput("awk '{print $5}' tests.txt").splitlines()

	# write relevant tests into lists
	# note: tests in txt file are ordered from newest to oldest, so only add to list if not already in there
	relevant_ids = []
	relevant_name = []
	for i in range(0, len(testid)):
	
		# if test not already in list and is not debug test
		# OTHER FILTERING CRITERIA: ADD HERE
		if name[i] not in relevant_name and not name[i].endswith("debug"):
			if (Options.id != "None" and testid[i] == Options.id) or (Options.sb == "True" and "sb_" in name[i]) or (Options.sb == "False" and "sb_" not in name[i]):

				relevant_ids.append( testid[i])
				relevant_name.append( name[i] )

	# turn lists into list of dicts
	tests = []
	for j in range(0, len(testid)):

		# if test is relevant, add to tests list
		if testid[j] in relevant_ids:
			tmp = {}
			tmp["test_id"] = testid[j]
			tmp["name"] = name[j]
			tmp["status"] = teststatus[j]
			tmp["revision"] = revision[j]
			tmp["dir"] = testdir[j]
			tests.append( tmp )
	
	return tests

################################################
def download_test_files(prefix):

#	s = get_url('/summaries/master').json()

#	for k in s: print(k)
#	print(f'Summary for branch master:\n', json.dumps(s, sort_keys=True, indent=2))
#	print(f'Tests for branch master:\n', json.dumps(s['tests'], sort_keys=True, indent=2))

	# get tests to download, this is a list of dicts:
	# [{'test_id': '674827', 'name': 'scientific.sb_talaris14_relax_fast_5iter', 'status': 'failed', 'dir': 'https://b3.graylab.jhu.edu/test/674827'}, {'test_id': '674826', 'name': 'scientific.sb_talaris14_relax_fast', 'status': 'failed', 'dir': 'https://b3.graylab.jhu.edu/test/674826'}]
	# put filtering criteria into the get_test_ids() function
	tests = get_test_ids()
	
#	tests = [{'test_id': '674827', 'name': 'scientific.sb_talaris14_relax_fast_5iter', 'status': 'failed', 'dir': 'https://b3.graylab.jhu.edu/test/674827'}]

	scientific_tests = [
		t for t in tests
#		if t['name'].startswith('scientific.') and not t['name'].endswith('.debug') and not t['name'].startswith('scientific.protein_data_bank_diagnostic')
#		if t['name'].startswith('scientific.') and not t['name'].endswith('.debug')
	]
	
	scientific_test_names = [ t['name'] for t in scientific_tests ]
	print(f'Scientific tests: ', scientific_test_names)

	url = urllib.parse.urlparse(_server_url_)
	user, password = Config.get('main', 'user'), Config.get('main', 'password')
	port = f':{url.port}' if url.port else ''

	all_tests = []

	logfile = 'log_' + datetime.today().strftime('%Y-%m-%d')
	print (logfile)
	with open( logfile, 'w' ) as f:

		for test in scientific_tests:
			
			# if test id given, only get that test
			# if test id not given, this if statement won't hold and all tests will be downloaded
			if Options.id is not None and test["test_id"] != Options.id:
				continue
			
			print ("downloading data for", test["name"], "with", test["test_id"])
			
			# create directories on local machine
			path = f'{prefix}/{test["name"]}'
			if not os.path.exists(path): os.makedirs(path)

			# create command line for download
#			command_line = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-host-directories --cut-dirs=5 {url.scheme}://{url.hostname}{port}/api/v1/tests/{test["test_id"]}/files/index.html'
#			command_line = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}{port}/api/v1/tests/{test["test_id"]}/files'
#			command_line = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}/test/{test["test_id"]}/file-list'
			command_line = f'wget \'{url.scheme}://{url.hostname}/test/{test["test_id"]}/file-list\' --user {user} --password {password} --no-check-certificate -O test_files'

			print (command_line)

			# download json file with dictionary containing list of files to download
			os.system( command_line )
			cmd = 'grep href test_files | awk -F\\" \'{print $2}\''
			files = subprocess.getoutput( cmd ).splitlines()

			# read file list
			for fn in files:
				
				# download each of those files
			
				if "benchmark" in fn:
					continue
				
				if "log" in fn:
					continue

				if "data" in fn:
					continue
				
#				if fn.endswith(".pdb"):
#					continue

#				if fn.endswith(".pdb.gz"):
#					continue

				if fn.endswith(".in_progress"):
					continue

				if fn.endswith(".condor"):
					continue

				if fn.endswith("index.html") or fn.endswith(".png"):
					
					print ("===", fn, flush=True)
				
#				cmd = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}/' + fn

				if Options.pdbid.upper() in fn:
					cmd = f'cd {path} && wget --user {user} --password {password} --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}/' + fn
					print (cmd)
					os.system( cmd )

				# to PDF files, use weasyprint outside of python
				# > weasyprint index.html out.pdf
#				if fn.endswith("index.html"):
#					print (file)
#					summary = f'{path}/' + f.split("/")[-1]
#					outfile = 'log_' + test["name"] + '.pdf'
#					os.system( "weasyprint index.html " + outfile )
#					pdfkit.from_file(summary, outfile)
#					HTML( summary ).write_pdf( outfile )

			# write revisions file
			with open( path + '/revision', 'w' ) as fr:
				if ("revision" in test):
					fr.write( test["revision"] + "\n" )
				fr.write( "test_id:" + test["test_id"] + "\n" )
			fr.close()
			
			# add revision and test number to end of readme
			os.system( "echo '\n## REVISION' >> " + path + "/readme.md" )
			os.system( "cat " + path + "/revision >> " + path + "/readme.md" )		

	f.close()

################################################
################################################

def main(args) -> None:
	''' scientific_test_report.py main '''

	parser = ArgumentParser(description=main.__doc__)

	parser.add_argument("--config", default=None, action="store", help="Location of .ini file with main configuration. Default is to use `config.ini` from the script dir.")
	parser.add_argument("--id", default=None, action="store", help="Test id to get files for. This is the last number in the URL https://b3.graylab.jhu.edu/test/555956 that you can get from the Benchmark server.")
	parser.add_argument("--pdbid", default=None, action="store", help="PDBID to get files for.")
	parser.add_argument("--sb", default=False, action="store", help="Download scorefunction comparison tests, yes or no. They start with sb_, therefore the flag.")
	
	global Options
	Options = parser.parse_args()
	global test_number
	test_number = Options.id

	if Options.id is None:
#		sys.exit( "ERROR: test id not given - please provide one with '--id <test_id>'. Exiting." )
		print ("WARNING: No test id given. Using the latest tests in master.")
		print ("NOTE: This is a lot of data!!!")
		test_number = "0000000"

	if Options.config is None: Options.config = os.path.dirname( os.path.abspath(__file__) ) + '/config.ini'

	global Config
	Config = ConfigParser( dict(here=os.path.abspath( os.path.dirname(Options.config ) ) ) )
	with open(Options.config) as f: Config.read_file(f)

	global _server_url_
	_server_url_ = Config.get('main', 'server.url')

	scientific_test_report_prefix = os.path.abspath( Config.get('DEFAULT', 'scientific_test_report_prefix') ) + "/" + datetime.today().strftime('%Y-%m-%d') + "_" + test_number
	 
#	if os.path.isdir(scientific_test_report_prefix): shutil.rmtree(scientific_test_report_prefix)  # always cleanup report dir
#	os.makedirs(scientific_test_report_prefix)

	download_test_files(scientific_test_report_prefix)


if __name__ == "__main__":
	main(sys.argv)
