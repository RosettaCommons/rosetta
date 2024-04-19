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

    ##################################################################################
    # Long list of double-checks to make sure we don't mess up the documentation state.

    branches = find_branches(rosetta_dir)
    if len(branches) != 1:
        return FAIL(None, f"Could not figure out which branch to update documentation for, commit belongs to following (expecting exactly one): " + ', '.join(branches) )
    main_branch = branches.pop()

    # Now update to the branch tip
    res, output = execute('Checking out branch...', 'cd {} && git fetch && git update-ref refs/heads/{branch} origin/{branch} && git reset --hard {branch} && git checkout {branch} && git submodule update && git branch --set-upstream-to=origin/{branch} {branch}'.format(rosetta_dir, branch=main_branch), return_='tuple')
    if res:
        return FAIL(output, "Can't update to proper branch value.")


    res, output = execute("Getting main's revision number ...", f'cd {rosetta_dir} && git rev-parse HEAD', return_='tuple')
    if res:
        return FAIL(output, "Can't get main's revision number")
    main_rev = output.strip()

    res, output = execute("Getting documentation's revision number ...", f'cd {doc_dir} && git rev-parse HEAD', return_='tuple' )
    if res:
        return FAIL(output, "Can't get documentation's revision number")
    doc_rev = output.strip()

    # I think this should probably work decently.
    doc_branch = main_branch

    res, output = execute("Does the doc branch exist?", f'cd {doc_dir} && git rev-parse --verify origin/{doc_branch}', return_='tuple')
    if res:
        # Branch does not exist - create it.
        res, output = execute("Creating branch for documentation updates", f'cd {doc_dir} && git checkout -b {doc_branch} && git branch --set-upstream-to=origin/{doc_branch} {doc_branch}', return_='tuple' )
    else:
        # Branch already exists - check it out, and then make sure we include what we have currently for the documentation
        res, output = execute("Checking out branch for documentation updates", f'cd {doc_dir} && git checkout {doc_branch} && git pull && git merge {doc_rev}', return_='tuple')
        if res:
            return FAIL(output, f"Can't properly checkout documentation repo branch {doc_branch}")

    ##########################################################
    # Update the repo(s)

    # try: finally: here to make sure we clean up the documentation directory in case of failures
    try:

        # Compile Rosetta - will also update the options documentation.
        res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)

        if res:
            results[_StateKey_] = _S_build_failed_
            results[_LogKey_]   = f'Compiling: {build_command_line}\n{output}'
            return results

        # Now attempt to generate the XSD

        ext = calculate_extension(platform, mode)
        res, output =  execute('Generating RosettaScripts XSD ...', f'cd {doc_dir} && {rosetta_dir}/source/bin/rosetta_scripts.{ext} -parser::output_schema rosettascripts.xsd', return_='tuple' )
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue running XSD generation: \n' + output
            return results


        # Clear out the old md files in the XSD directory
        for f in glob.glob(f"{doc_dir}/scripting_documentation/RosettaScripts/xsd/*.md"):
            os.remove(f)

        # Parse the XSD into the appropriate format.
        res, gen_output = execute('Generating Markdown from XSD ...', f'cd {doc_dir} && python2 {rosetta_dir}/tools/doc_tools/xsd_to_doc_fragments.py rosettascripts.xsd  ./scripting_documentation/RosettaScripts/xsd/', return_='tuple')
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue converting XSD to Markdown: \n' + gen_output
            return results

        # We need to double check if there's any case sensitivity issues
        res, output = execute('Checking for case sensitivity issues ...', f'cd {doc_dir} && find . | sort -f | uniq -di', return_='tuple')
        if res or output.strip() != '':
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Case sensitivity failure encountered: \n' + gen_output + '\n\n' + output
            return results

        # While we're at it, we can update the full options list - should be regenerated from the compile
        shutil.copyfile(f'{rosetta_dir}/source/src/basic/options/full-options-list.md', f'{doc_dir}/full-options-list.md')

        #################################
        # Cleanup

        os.remove(f"{doc_dir}/rosettascripts.xsd")

        ####################
        # Commit to repo

        # Add XSD changes to staging area
        execute('Adding Markdown changes to XSD ...', f'cd {doc_dir} && git add ./scripting_documentation/RosettaScripts/xsd/' )
        execute('Removing any removed Markdown XSD files...', f'cd {doc_dir} && git add -u ./scripting_documentation/RosettaScripts/xsd/' )

        # Don't bother adding full option list unless more than the date changes
        res, output = execute('Checking fulloption list differences ...', f'cd {doc_dir} && git --no-pager diff --no-color --numstat -- full-options-list.md', return_='tuple')
        if len(output.split()) == 3 and output.split() == ['1','1','full-options-list.md']:
            execute('Reverting date-only change to full-options-list.md ...', f'cd {doc_dir} && git checkout -- ./full-options-list.md' )
        else:
            # Add that to the staging area
            execute('Add Option list update ...', f'cd {doc_dir} && git add ./full-options-list.md' )

        res, status_output = execute('Getting status output ...', f'cd {doc_dir} && git status', return_='tuple')

        # Now we commit
        res, output = execute('Committing ...', f'cd {doc_dir} && git commit -m "Test server update for {main_rev}"', return_='tuple' )
        if res:
            # We get this when there's nothing that needs to be updated.
            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = 'No update to commit: \n' + output + '\n'
            return results

        # Push to Github documentation repo
        res, output = execute('Pushing ...', f'cd {doc_dir} && git push', return_='tuple'  )
        if res:
            return FAIL(output, "Couldn't push documentation update")

        # Update main and push to branch on main
        res, output = execute("Updating current branch", f'cd {rosetta_dir} && git add documentation/ && git commit -m "Test server updating documentation submodule"', return_='tuple')
        if res:
            return FAIL(output, "Couldn't update main repo")
        res, output = execute("Pushing main...", f'cd {rosetta_dir} && git push', return_='tuple' )
        if res:
            return FAIL(output, "Couldn't push main repo updates")

        # If we've reached here, we're successful
        results[_StateKey_] = _S_passed_
        results[_LogKey_]   = 'Updated documentation repo : \n' + output + '\n\n' + gen_output + '\n\n' + status_output
        return results

    finally:
        # Revert any changes that haven't been properly pushed, also remove any extra (non-ignored) files
        res, output = execute("What should the doc repo be at?...", f'cd {rosetta_dir} && git ls-tree HEAD -- documentation/', return_='tuple' )
        reset_to = output.split()[2]
        execute('Cleaning up', f'cd {doc_dir} && ( git reset --hard {reset_to}; git clean -fd )' )

def run(test, repository_root, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if test == 'documentation': return run_documentation_update(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknown maintenace test: {}!'.format(test))
