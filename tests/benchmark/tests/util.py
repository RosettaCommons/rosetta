#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/util.py
## @brief  Benchmark utility tests for various auxiliary tasks
## @author Sergey Lyskov

import os, os.path, json, shutil, stat, glob, collections

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'  # api version

def generate_rosetta_app_list(rosetta_dir, working_dir, platform, config, hpc_driver, verbose, debug):

    scons_app_files = { 'apps.src.settings' : 'public', 'pilot_apps.src.settings.all' : 'pilot' }

    apps = collections.defaultdict(list)

    for scons_app_file, tag in scons_app_files.items():
        with open(rosetta_dir + '/source/src/' + scons_app_file) as f:
            code = compile(f.read(), scons_app_file, 'exec')
            globals = {}
            exec(code, globals)

            for dir_ in globals['sources']:
                for f in globals['sources'][dir_]:
                    apps[tag].append( f.split('/')[-1] )

    for v in apps.values(): v.sort()

    results = dict( apps = apps, revision = dict( branch = config['branch'], revision_id = config['revision'] ) )

    with open( working_dir + '/rosetta-apps.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)

    return {_StateKey_  : _S_passed_, _ResultsKey_ : results, _LogKey_ : '' }


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test =='apps': return generate_rosetta_app_list(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError(f'Unknown scripts test: {test}!')
