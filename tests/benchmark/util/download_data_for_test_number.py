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
import pdfkit
from datetime import datetime


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

def get_url (name, **args):   return http_request(requests.get,	name, **args)
def post_url(name, **args):   return http_request(requests.post,   name, **args)
def delete_url(name, **args): return http_request(requests.delete, name, **args)



def download_test_files(prefix):
	s = get_url('/summaries/master').json()

#	print (s)

	tests = [
		{"test_id": "599753", "name": "scientific.sb_talaris14_relax_fast_5iter"},
		{"test_id": "599752", "name": "scientific.sb_talaris14_relax_fast"},
		{"test_id": "599751", "name": "scientific.sb_talaris14_relax_cartesian"},
		{"test_id": "599750", "name": "scientific.sb_talaris14_loop_modeling_ngk_12res"},
		{"test_id": "599749", "name": "scientific.sb_talaris14_loop_modeling_kic_fragments_12res"},
		{"test_id": "599748", "name": "scientific.sb_talaris14_loop_modeling_kic_12res"},
		{"test_id": "599747", "name": "scientific.sb_talaris14_loop_modeling_ccd_12res"},
		{"test_id": "599746", "name": "scientific.sb_talaris14_fast_design"},
		{"test_id": "599745", "name": "scientific.sb_talaris14_docking"},
		{"test_id": "599744", "name": "scientific.sb_talaris13_relax_fast_5iter"},
		{"test_id": "599743", "name": "scientific.sb_talaris13_relax_fast"},
		{"test_id": "599742", "name": "scientific.sb_talaris13_relax_cartesian"},
		{"test_id": "599741", "name": "scientific.sb_talaris13_loop_modeling_ngk_12res"},
		{"test_id": "599740", "name": "scientific.sb_talaris13_loop_modeling_kic_fragments_12res"},
		{"test_id": "599739", "name": "scientific.sb_talaris13_loop_modeling_kic_12res"},
		{"test_id": "599738", "name": "scientific.sb_talaris13_loop_modeling_ccd_12res"},
		{"test_id": "599737", "name": "scientific.sb_talaris13_fast_design"},
		{"test_id": "599736", "name": "scientific.sb_talaris13_docking"},
		{"test_id": "599735", "name": "scientific.sb_score12_relax_fast_5iter"},
		{"test_id": "599734", "name": "scientific.sb_score12_relax_fast"},
		{"test_id": "599733", "name": "scientific.sb_score12_relax_cartesian"},
		{"test_id": "599732", "name": "scientific.sb_score12_loop_modeling_ngk_12res"},
		{"test_id": "599731", "name": "scientific.sb_score12_loop_modeling_kic_fragments_12res"},
		{"test_id": "599730", "name": "scientific.sb_score12_loop_modeling_kic_12res"},
		{"test_id": "599729", "name": "scientific.sb_score12_loop_modeling_ccd_12res"},
		{"test_id": "599728", "name": "scientific.sb_score12_fast_design"},
		{"test_id": "599727", "name": "scientific.sb_score12_docking"},
		{"test_id": "599726", "name": "scientific.sb_ref2015_relax_fast_5iter"},
		{"test_id": "599725", "name": "scientific.sb_ref2015_relax_fast"},
		{"test_id": "599724", "name": "scientific.sb_ref2015_relax_cartesian"},
		{"test_id": "599723", "name": "scientific.sb_ref2015_loop_modeling_ngk_12res"},
		{"test_id": "599722", "name": "scientific.sb_ref2015_loop_modeling_kic_fragments_12res"},
		{"test_id": "599721", "name": "scientific.sb_ref2015_loop_modeling_kic_12res"},
		{"test_id": "599720", "name": "scientific.sb_ref2015_loop_modeling_ccd_12res"},
		{"test_id": "599719", "name": "scientific.sb_ref2015_fast_design"},
		{"test_id": "599718", "name": "scientific.sb_ref2015_docking"},
		{"test_id": "599717", "name": "scientific.sb_ligand_docking"}
	]

	tests = [
		{"test_id": "599728", "name": "scientific.sb_score12_fast_design"},
		{"test_id": "605064", "name": "scientific.sb_score12_loop_modeling_kic_fragments_12res"},
		{"test_id": "604747", "name": "scientific.sb_ref2015_loop_modeling_kic_fragments_12res"},
		{"test_id": "604748", "name": "scientific.sb_ref2015_loop_modeling_ngk_12res"},
	]

	tests = [
		{"test_id": "619450", "name": "scientific.sb_score12_fast_design"}
	]
		
#	for k in s: print(k)
#	print(f'Summary for branch master:\n', json.dumps(s, sort_keys=True, indent=2))
#	print(f'Tests for branch master:\n', json.dumps(s['tests'], sort_keys=True, indent=2))

	scientific_tests = [
		t for t in tests
#		if t['name'].startswith('scientific.') and not t['name'].endswith('.debug') and not t['name'].startswith('scientific.protein_data_bank_diagnostic')
#		if t['name'].startswith('scientific.') and not t['name'].endswith('.debug')
	]

	print (scientific_tests)

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
			
#			if test["test_id"] != int(Options.id):
#				continue

			print ("downloading data for", test["name"], "with", test["test_id"])
			
			path = f'{prefix}/{test["name"]}'
			if not os.path.exists(path): os.makedirs(path)

			# create command line for download
#			command_line = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-host-directories --cut-dirs=5 {url.scheme}://{url.hostname}{port}/api/v1/tests/{test["test_id"]}/files/index.html'
			command_line = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}{port}/api/v1/tests/{test["test_id"]}/files'

			# download json file with dictionary containing list of files to download
			output = subprocess.getoutput(command_line)
			print (output)

			# read json file
			with open( path + "/files", 'r' ) as json_file:
				test_info = json.load( json_file )
				
				# download each of those files
				for file in test_info['files']:
					
					if file.startswith("benchmark"):
						continue
					
					if file.startswith("hpc-logs"):
						continue

					if file.endswith(".pdb"):
						continue

					if file.endswith(".pdb.gz"):
						continue

					if file.endswith(".in_progress"):
						continue

					if file.endswith(".condor"):
						continue
					
					cmd = f'cd {path} && wget --user {user} --password {password} --recursive --directory-prefix={path} --no-parent --no-check-certificate --cut-dirs=5 -nH {url.scheme}://{url.hostname}{port}/api/v1/tests/{test["test_id"]}/files/' + file
					
					print (cmd)
					
					output = subprocess.getoutput(cmd)

	f.close()


def main(args) -> None:
	''' scientific_test_report.py main '''

	parser = ArgumentParser(description=main.__doc__)

	parser.add_argument("--config", default=None, action="store", help="Location of .ini file with main configuration. Default is to use `config.ini` from the script dir.")
	parser.add_argument("--id", default=None, action="store", help="Test id to get files for. This is the last number in the URL https://b3.graylab.jhu.edu/test/555956 that you can get from the Benchmark server.")
	
	global Options
	Options = parser.parse_args()

	if Options.id is None:
		sys.exit( "ERROR: test id not given - please provide one with '--id <test_id>'. Exiting." )

	if Options.config is None: Options.config = os.path.dirname( os.path.abspath(__file__) ) + '/config.ini'

	global Config
	Config = ConfigParser( dict(here=os.path.abspath( os.path.dirname(Options.config ) ) ) )
	with open(Options.config) as f: Config.read_file(f)

	global _server_url_
	_server_url_ = Config.get('main', 'server.url')

	scientific_test_report_prefix = os.path.abspath( Config.get('DEFAULT', 'scientific_test_report_prefix') ) + "/" + datetime.today().strftime('%Y-%m-%d') + "_" + Options.id
	 
#	if os.path.isdir(scientific_test_report_prefix): shutil.rmtree(scientific_test_report_prefix)  # always cleanup report dir
#	os.makedirs(scientific_test_report_prefix)

	download_test_files(scientific_test_report_prefix)


if __name__ == "__main__":
	main(sys.argv)
