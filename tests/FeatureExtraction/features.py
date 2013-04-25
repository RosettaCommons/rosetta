# -*- coding: utf-8 -*-
# :noTabs=true:

#
# Full log-output of this scripts will be saved by FeatureExtraction daemon to .features.py.log file for test FeatureExtraction for corresponding revision
#

import os.path

try: import json
except ImportError: import simplejson as json

args = json.load( file('.arguments.py') )

print 'Features script, running with args:', json.dumps(args, sort_keys=True, indent=2)
working_dir = args['working_dir']


# list of tags telling us why this test was scheduled to run, it could contain following: 'Plain', 'ScoreFunctionChange', 'Release'
# this field will be automatically added by damon to output.json
tags = args['tags']

# dict containing list of FeaturesDB revisions avaliable for reference. Each revision is just plain integer
# files for each revisions could be accessed from: args['reference_DB']/<rev_number>
reference_info = args['reference_info']

reference_DB = args['reference_DB']
# use reference_DB/<rev_number> to obtain information from previous runs. Each dir preserved as-it was by the end of this script run.
# be carreful to not modify contents of reference_DB/... - it is not guarantee to be 'fresh' for each run. You may however unpack SQLite DB if needed.

# Writing some sample features results
results_ExtractionDB = args['results_ExtractionDB']
f = file(results_ExtractionDB + 'readme.txt', 'w');  f.write('Write feature DB results files here!');  f.close()


to_delete = []  # list of revisions that should be removed from Features DB after this finish running. Only revision from args['reference_info'] could be used

# Writing output.json, it will be stored by FeaturesExtraction daemon in to FeatureExtraction.yaml field
# root element _must_ be a dictionary. Following fields is reserved and should not be used: 'tags', 'state'
output = dict(to_delete=to_delete)

json.dump(output, file(results_ExtractionDB + 'output.json', 'w') )




# Writing some sample report results
results_FeaturesPlots = args['results_FeaturesPlots']
f = file(results_FeaturesPlots + 'index.html', 'w');  f.write('<html><body>Sample Report!</body></html>');  f.close()

# Writing output.json, for plots. it will be stored by FeaturesExtraction daemon in to Features.yaml field
# content of it will be used to determent of results avalability
output = dict()

json.dump(output, file(results_FeaturesPlots + 'output.json', 'w') )

