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

import os, os.path, json, commands, shutil, stat, glob

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'  # api version

def run_documentation_update(rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    memory = config['memory']
    jobs = config['cpu_count']
    mode = 'release'
    doc_dir = os.path.join( rosetta_dir, '../documentation/' )
    results = {}

    TR = Tracer(verbose)

    TR('Running Documentation Maintenance "test" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    # Check to make sure that we're appropriately set up
    if config['branch'] not in 'master commits release unknown':
        results[_StateKey_] = _S_script_failed_
        results[_LogKey_]   = "Can't run documentation maintenance with non-master testing branch. Found testing branch: {}\n".format(config['branch'])
        return results

    res, output = execute("Getting main's revision number ...", 'cd {} && git rev-parse HEAD'.format(rosetta_dir), return_='tuple')
    if res:
        results[_StateKey_] = _S_script_failed_
        results[_LogKey_]   = "Can't get main's revision numbers: {}\n".format(output.strip())
        return results
    main_rev = output.strip()

    res, output = execute('Checking documentation branch name ...', 'cd {} && git rev-parse --abbrev-ref HEAD'.format(doc_dir), return_='tuple' )
    doc_branch_name = output.strip()
    if res or not ( doc_branch_name == 'master' or doc_branch_name == 'roccomoretti/benchmark_doc_gen'):
        results[_StateKey_] = _S_script_failed_
        results[_LogKey_]   = "Can't run documentation maintenance with documentation not on branch master. Found branch: {}\n".format(output.strip())
        return results

    res, output = execute('Checking that documentation repo is clean ...', 'cd {} && git status --porcelain'.format(doc_dir), return_='tuple' )
    if res or output.strip() != '':
        results[_StateKey_] = _S_script_failed_
        results[_LogKey_]   = "Can't run documentation maintenance with unclean documentation repo: \n\n" + output
        return results

    # We assume that the test system hooks have already updated both branches to the latest from GitHub.

    # Compile Rosetta - will also update the options documentation.
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, config, mode=mode, verbose=verbose)

    if res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Compiling: {}\n'.format(build_command_line) + output
        return results

    # try: finally: here to make sure we clean up the documentation directory in case of failures

    try:
        # Now attempt to generate the XSD

        ext = calculate_extension(platform, mode)
        res, output =  execute('Generating RosettaScripts XSD ...', 'cd {doc_dir} && {rosetta_dir}/source/bin/rosetta_scripts.{ext} -parser::output_schema rosettascripts.xsd'.format(**vars()), return_='tuple' )
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue running XSD generation: \n' + output
            return results


        # Clear out the old md files in the XSD directory
        for f in glob.glob("{doc_dir}/scripting_documentation/RosettaScripts/xsd/*.md".format(**vars())):
            os.remove(f)

        # Parse the XSD into the appropriate format.
        res, gen_output = execute('Generating Markdown from XSD ...', 'cd {doc_dir} && {rosetta_dir}/../tools/doc_tools/xsd_to_doc_fragments.py rosettascripts.xsd  ./scripting_documentation/RosettaScripts/xsd/'.format(**vars()), return_='tuple')
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue converting XSD to Markdown: \n' + gen_output
            return results

        # We need to double check if there's any case sensitivity issues
        res, output = execute('Checking for case sensitivity issues ...', 'cd {doc_dir} && find . | sort -f | uniq -di'.format(**vars()), return_='tuple')
        if res or output.strip() != '':
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Case sensitivity failure encountered: \n' + gen_output + '\n\n' + output
            return results

        # Now we can add to staging area
        execute('Adding Markdown changes to XSD ...', 'cd {doc_dir} && git add ./scripting_documentation/RosettaScripts/xsd/'.format(**vars()) )
        execute('Removing any removed Markdown XSD files...', 'cd {doc_dir} && git add -u ./scripting_documentation/RosettaScripts/xsd/'.format(**vars()) )

        os.remove("{doc_dir}/rosettascripts.xsd".format(**vars()))

        # While we're at it, we can update the full options list - should be regenerated from the compile
        shutil.copyfile('{rosetta_dir}/source/src/basic/options/full-options-list.md'.format(**vars()), '{doc_dir}/full-options-list.md'.format(**vars()))

        # Don't bother adding full option list unless more than the date changes
        res, output = execute('Checking fulloption list differences ...', 'cd {doc_dir} && git diff --numstat -- full-options-list.md'.format(**vars()), return_='tuple')
        if len(output.split()) == 3 and output.split() == ['1','1','full-options-list.md']:
            execute('Reverting date-only change to full-options-list.md ...', 'cd {doc_dir} && git checkout -- ./full-options-list.md'.format(**vars()) )
        else:
            # Add that to the staging area
            execute('Add Option list update ...', 'cd {doc_dir} && git add ./full-options-list.md'.format(**vars()) )

        res, status_output = execute('Getting status output ...', 'cd {doc_dir} && git status'.format(**vars()), return_='tuple')

        # Now we commit
        res, output = execute('Committing ...', 'cd {doc_dir} && git commit -m "Test server update for {main_rev}"'.format(**vars()), return_='tuple' )
        if res:
            # We get this when there's nothing that needs to be updated.
            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = 'No update to commit: \n' + output + '\n'
            return results

        # Push to master
        res, output = execute('Pushing ...', 'cd {doc_dir} && git push'.format(**vars()).format(**vars()), return_='tuple'  )
        if res == 0:
            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = 'Updated documentation repo : \n' + output + '\n\n' + gen_output + '\n\n' + status_output
            return results

        # Merge issues? Try merging and retrying.
        res, output = execute('Pulling ...', 'cd {doc_dir} && git pull --no-edit'.format(**vars()).format(**vars()), return_='tuple'  )
        if res:
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue merging with GitHub : \n' + output
            return results

        # And re-push
        res, output = execute('Re-Pushing ...', 'cd {doc_dir} && git push'.format(**vars()).format(**vars()), return_='tuple'  )
        if res:
            # If the push fails twice, there's a good chance something else is going on.
            results[_StateKey_] = _S_script_failed_
            results[_LogKey_]   = 'Issue pushing to GitHub : \n' + output
        else:
            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = 'Updated documentation repo : \n' + output + '\n\n' + gen_output + '\n\n' + status_output
            return results

    finally:
        # Revert any changes that haven't been properly pushed, also remove any extra (non-ignored) files
        # Remember that we checked that we were on the documentation directory's master branch
        execute('Cleaning up', 'cd {doc_dir} && ( git reset --hard origin/{doc_branch_name}; git clean -fd )'.format(**vars()) )

def run(test, rosetta_dir, working_dir, platform, config, hpc_driver, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if test == 'documentation': return run_documentation_update(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknown maintenace test: {}!'.format(test))

