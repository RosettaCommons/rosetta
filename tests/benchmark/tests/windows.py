#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/windows.py
## @brief  Benchmark tests for various tests for Windows platform
## @author Sergey Lyskov

import os, os.path, json, shutil, stat, glob, collections, codecs

from subprocess import getstatusoutput, getoutput

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'  # api version

#_activate_cl_ = 'C:\\BuildTools\\Common7\\Tools\\VsDevCmd.bat -arch=x64 -host_arch=amd64'
_activate_cl_ = '%ACTIVATE_MSVC%'

def run_git_clean_if_needed(rosetta_dir):
    signature_file = rosetta_dir + '/.windows-vc-tools-signature'

    expected_signature = getoutput(f'{_activate_cl_} && set VCToolsVersion && cl')

    if os.path.isfile(signature_file):
        with open(signature_file) as f: signature = f.read()
    else: signature = ''

    if signature != expected_signature:
        res, output = getstatusoutput( f'cd {rosetta_dir} && git clean -dfx source/cmake' )
        print(output)
        if res: print('Git clean failed, aborting...'); sys.exit(1)

        with open(signature_file, 'w') as f: f.write(expected_signature)



def build_rosetta(mode, rosetta_dir, working_dir, platform, config, hpc_driver, verbose, debug):
    cpu_count = config['cpu_count']

    run_git_clean_if_needed(rosetta_dir)

    mode = dict(release='windows', debug='windows_debug')[mode]

    command_line = f'{_activate_cl_} && cd {rosetta_dir}\\source && python ninja_build.py {mode} -k -remake -j'

    print('Compiling...', command_line)
    res, output = getstatusoutput( f'{command_line}{cpu_count}' )

    if res: res, output = getstatusoutput( f'{command_line}1' ) # re-running with a single build thread to make log easier to read

    output = f'Running: {command_line}{cpu_count}\n{output}\n'

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(output)

    res_code = _S_failed_ if res else _S_passed_

    return {_StateKey_  : res_code, _ResultsKey_ : {}, _LogKey_ : output}


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test =='build.debug':   return build_rosetta('debug', repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if test =='build.release': return build_rosetta('release', repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError(f'Unknown scripts test: {test}!')
