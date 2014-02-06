#!/usr/bin/env python

import sys, argparse
import helpers, schema

# This script uses two sets of global variables.  The first are the command 
# line arguments, which are defined just below.  The second are the data 
# arrays, which are defined in the build_arrays() method.

parser = argparse.ArgumentParser()
helpers.add_job_arguments(parser)
parser.add_argument('--force', '-f', action='store_true')
arguments = parser.parse_args()

if __name__ == '__main__':
    job = arguments.job
    url = arguments.database

    with helpers.connect_to_database(url):
        query = schema.numpy_cache.query

        all_arrays_exist = (
                query.filter_by(job_id=job, type='iterations').count() and
                query.filter_by(job_id=job, type='scores').count() and
                query.filter_by(job_id=job, type='rmsds').count() and
                query.filter_by(job_id=job, type='torsions').count())

        if all_arrays_exist and not arguments.force:
            print "No work to do."
            sys.exit()
        
        helpers.cache_trajectory(job, arguments.force)

