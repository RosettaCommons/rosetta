#!/usr/bin/env python

import argparse
import helpers

parser = helpers.define_args()
arguments = parser.parse_args()

with helpers.connect_to_database(arguments.database):
    helpers.wipe_cache(arguments.job)

helpers.finish_script()
