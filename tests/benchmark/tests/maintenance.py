#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/maintenance.py
## @brief  Benchmark tests which are intended for housekeeping purposes
## @author Rocco Moretti (rmorettiase@gmail.com)

import os, os.path, json, shutil, stat, glob

# A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location
import importlib.util, sys
importlib.util.spec_from_file_location(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py').loader.exec_module(sys.modules[__name__])

_api_version_ = '1.1'

def FAIL(output, message):
    results = {}
    results[_StateKey_] = _S_script_failed_
    if output is None:
        results[_LogKey_]   = message
    else:
        output = output.strip()
        results[_LogKey_]   = f"{message}:\n{output}\n"
    return results


def find_branches(rosetta_dir):
    branches = set()
    o = execute('Getting list of branches where HEAD is present...', f'cd {rosetta_dir} && git branch --no-color -r --contains HEAD', return_='output')
    for line in o.split('\n'):
        line = line.replace('*', '').split()
        if not line:
            continue
        br = line[0]
        if br.startswith('origin/'):
            branches.add( br[len('origin/'):] )

    if 'main' in branches: branches = set( ['main'] )
    return branches

def run_documentation_update(rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    memory = config['memory']
    jobs = config['cpu_count']
    mode = 'release'
    doc_dir = os.path.join( rosetta_dir, 'documentation/' )
    results = {}

    TR = Tracer(verbose)

    TR('Running Documentation Maintenance "test" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    #### We assume that the documentation submodule has been checked out and updated to the current state, and the running directory is "clean"

    # try: finally: here to make sure we clean up the documentation directory afterward.
    try:

        # Compile Rosetta - will also update the options documentation.
        res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)

        if res:
            results[_StateKey_] = _S_build_failed_
            results[_LogKey_]   = f'Compiling: {build_command_line}\n{output}'
            return results

        # Update the options documentation -- should be regenerated from the compile.
        shutil.copyfile(f'{rosetta_dir}/source/src/basic/options/full-options-list.md', f'{doc_dir}/full-options-list.md')

        # Now attempt to generate the XSD

        ext = calculate_extension(platform, mode)
        res, xsd_output =  execute('Generating RosettaScripts XSD ...', f'cd {doc_dir} && {rosetta_dir}/source/bin/rosetta_scripts.{ext} -parser::output_schema rosettascripts.xsd', return_='tuple' )
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue running XSD generation: \n' + xsd_output
            return results

        # Clear out the old md files in the XSD directory
        for f in glob.glob(f"{doc_dir}/scripting_documentation/RosettaScripts/xsd/*.md"):
            os.remove(f)

        # Parse the XSD into the appropriate format.
        res, gen_output = execute('Generating Markdown from XSD ...', f'cd {doc_dir} && python {rosetta_dir}/tools/doc_tools/xsd_to_doc_fragments.py rosettascripts.xsd  ./scripting_documentation/RosettaScripts/xsd/', return_='tuple')
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue converting XSD to Markdown: \n' + gen_output
            return results

        # We need to double check if there's any case sensitivity issues
        res, case_output = execute('Checking for case sensitivity issues ...', f'cd {doc_dir} && find . | sort -f | uniq -di', return_='tuple')
        if res or case_output.strip() != '':
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Case sensitivity failure encountered: \n' + gen_output + '\n\n' + case_output
            return results

        ####################
        # Generate diff

        # Don't bother with full option list unless more than the date changes
        res, output = execute('Checking fulloption list differences ...', f'cd {doc_dir} && git --no-pager diff --no-color --numstat -- full-options-list.md', return_='tuple')
        if len(output.split()) == 3 and output.split() == ['1','1','full-options-list.md']:
            execute('Reverting date-only change to full-options-list.md ...', f'cd {doc_dir} && git checkout -- ./full-options-list.md' )

        # Cleanup unneeded file
        if os.path.exists(f"{doc_dir}/rosettascripts.xsd"):
            os.remove(f"{doc_dir}/rosettascripts.xsd")

        # Make sure that we pull in any new files from the XSD
        execute('Adding Markdown changes to XSD ...', f'cd {doc_dir} && git add ./scripting_documentation/RosettaScripts/xsd/' )

        res, porc_status_output = execute('Getting status output ...', f'cd {doc_dir} && git status --porcelain', return_='tuple')

        if len( porc_status_output.strip() ) == 0:
            # No changes -- report success.
            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = 'No documentation changes to report.'
            return results

        # Changes to report, create a patch file for users
        # (Can't automatically commit results to either repo anymore.)

        res, user_diff_output = execute('Getting status output ...', f'cd {doc_dir} && git --no-pager diff --no-color --stat=240 HEAD', return_='tuple')

        res_diff, diff = execute('Generating patch file for changes', f'cd {doc_dir} && git --no-pager diff --no-color HEAD', return_="tuple")
        with open(working_dir+'/full_diff.patch', 'w') as f: f.write(diff)

        output = '''Documentation Repository Needs Updating

See full_diff.patch in the `Test files` for the full list of changes.

You can apply these changes to your branch by downloading that patch file and executing

        cd rosetta/documentation/; git apply ~/Downloads/full_diff.patch

(You may need to use `dos2unix` on the full_diff.patch first.)

You will then need to commit and push those changes to the `documentation` repository (i.e. make a PR in the docs repo),
and then update your branch in the `rosetta` repo to point the documentation submodule to the updated commit.

Git status output of diff:
'''

        output += user_diff_output

        results[_StateKey_] = _S_failed_
        results[_LogKey_]   = output
        return results

    finally:
        # Revert any changes, also remove any extra (non-ignored) files
        execute('Cleaning up', f'cd {doc_dir} && ( git reset --hard HEAD; git clean -fd )' )

def run(test, repository_root, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if test == 'documentation': return run_documentation_update(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknown maintenace test: {}!'.format(test))
