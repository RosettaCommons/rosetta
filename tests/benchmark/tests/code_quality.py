#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   serialization.py
## @brief  scientific/serialization.py
## Benchmark script for running serialization test
## @author Sergey Lyskov

import os, json
import glob, shutil
import plistlib
import time

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from ../__init__ import *, but init path is calculated relatively to this location


_api_version_ = '1.0'  # api version


def run_serialization_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    jobs = config['cpu_count']

    execute('Updating options, ResidueTypes and version info...', 'cd {}/source && ./update_options.sh && ./update_ResidueType_enum_files.sh && python version.py'.format(rosetta_dir) )

    executable = install_llvm_tool('serialization_validator', source_location='{}/../tools/clang_ast_transform/rosetta-refactor-tool'.format(rosetta_dir), config=config)
    #executable_path = executable.rpartition('/')[0]

    tools_dir = os.path.abspath(rosetta_dir+'/../tools')

    json_output_dir = working_dir + '/serialization-validator-output'
    os.mkdir( json_output_dir )

    command_line = 'export PYTHONPATH={tools_dir}/python_cc_reader${{PYTHONPATH+:$PYTHONPATH}} && cd {rosetta_dir}/source && python ../../tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e '\
    '"python  {rosetta_dir}/../tools/clang_ast_transform/run_serialization_validator_on_file.py --executable_path {executable} --json_output_path {json_output_dir} --filename" -n {jobs}'.format(**vars())

    res, output = execute('Running...', command_line, return_='tuple')
    output = 'Running: ' + command_line + '\n' + output + '\n\n'

    json_files = [f for f in os.listdir(json_output_dir) if f.endswith('.json')]
    if json_files:
        state = _S_failed_ if res else _S_passed_
        tests = {}
        results = {_TestsKey_:tests}

        for f in json_files:
            with open(json_output_dir+'/'+f) as _fd: jr = json.load(_fd)
            for k in jr:
                key = k[ len('src/') : ].replace('/', '_')
                if _StateKey_ not in jr[k]  or  jr[k][_StateKey_] != _S_passed_:
                    tests[key] = jr[k]
                    state = _S_failed_
                else:
                    output += 'Test {} is passed!\n'.format(key)

    else:
        state, results = _S_script_failed_, {}
        output += 'ERROR: serialization_validator tool run but i could not find json files... Terminating with script failure!'

    return {_StateKey_ : state,  _ResultsKey_ : results,  _LogKey_ : output }



def parse_clangSA_exclusions( exclusion_list, verbose ):
    exclusions = {}

    if exclusion_list is not None and os.path.exists(exclusion_list):
        with open(exclusion_list) as f:
            for line in f:
                line = line.split(None,1) # Split on first whitespace separator
                if len(line) != 2 or line[0][0] == '#':
                    continue
                filename = line[0].rsplit(':',1)[0]
                # Since line numbers can change, we don't want to use them for exclusions
                message = line[1].strip()
                exclusions.setdefault(filename,[]).append(message)
    else:
        Tracer(verbose)("Exclusions file "+exclusion_list+" not found, ignoring.")

    return exclusions

def parse_clangSA_plist_file( filename, root_dir, exclusions, verbose ):
    issues = {}
    with open(filename,'rb') as f: # readPlist needs binary-formatted files
        try:
            plist = plistlib.readPlist(f)
        except:
            Tracer(verbose)("Issue parsing plist file: "+filename+"- ignoring.")
            return issues, None

    files = plist.get('files',[])
    clang_version = plist.get('clang_version',None)

    for diag in plist.get('diagnostics',[]):
        description = diag.get('description','NO DESCRIPTION').strip()
        html_files = diag.get('HTMLDiagnostics_files', [] )
        location = diag.get('location',{})
        line = location.get('line',0)
        filenum = location.get('file',-1)
        if 0 <= filenum < len(files):
            filename = os.path.realpath(files[filenum])
            filename = os.path.relpath(filename, os.path.realpath(root_dir))
            if filename.startswith("external") or filename.startswith("src/apps/pilot"):
                continue
            if filename in exclusions and description in exclusions[filename]:
                continue
            issues.setdefault((filename, line, description), []).extend(html_files)

    return issues, clang_version

def run_clang_analysis_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    jobs = config['cpu_count']
    start_time = time.time()
    TR = Tracer(verbose)
    full_log = ''
    results = {}

    TR('Running test_suite="{}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format(__name__, **vars() ) )

    run_command_line = 'cd {rosetta_dir}/source/cmake/build_clangSA/ && ./go.sh {jobs}'.format(**vars())
    run_res, run_output = execute('Running clang static analysis...', run_command_line, return_='tuple' )

    # We don't exit on clang static analysis failure here - that's deliberate
    # We want need to touch all the broken files we know about later on.

    full_log += run_output

    #Find the output directory
    out_dirs = list(glob.glob("{rosetta_dir}/source/build/clang_SA/*".format(**vars())))
    if len(out_dirs) < 1:
        # We use the --keep-empty option of scan-build, so a missing directory means something has gone wrong
        TR('No output directory found!')
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'No outdir when running clang static analysis...: {}\n\n'.format(run_command_line) + full_log
        return results

    outdir = sorted(out_dirs)[-1] # The most recent one
    try:
        # Copy the html output to the saved output directory
        for f in glob.glob("{}/*.html".format(outdir)):
            shutil.copy(f, working_dir)

        # Parse the plist output to find which files have issues
        exclusions = parse_clangSA_exclusions( "{}/tests/benchmark/util/clangSA_exclusions.txt".format(rosetta_dir), verbose )

        issues = {}
        unable_to_parse_plist = []
        clang_version = set()
        for filename in glob.glob('{}/*.plist'.format(outdir)):
            file_issues, file_version = parse_clangSA_plist_file(filename, os.path.join(rosetta_dir,'source'), exclusions, verbose)
            for k in file_issues:
                issues.setdefault(k,[]).extend(file_issues[k])
            if file_version is None:
                unable_to_parse_plist.append( filename )
            else:
                clang_version.add( file_version )
        clang_version = '//'.join( [v for v in clang_version if v is not None] )

    finally:
        shutil.rmtree( outdir )
        pass

    historical = parse_clangSA_exclusions( "{}/tests/benchmark/util/clangSA_historical.txt".format(rosetta_dir), verbose )

    json_results = dict(tests={}, summary=dict(total=[], failed=[], failed_tests=[]))
    problem_log = []
    historical_log = []
    for filename in unable_to_parse_plist:
        problem_log.append("Cannot parse plist file " + filename)
        json_results['summary']['failed_tests'].append(filename)
    for f, l, m in sorted(issues.keys()):
        # Update the modification date on files with problems, so they show up on the next run.
        # (If we do this multiple times per file, it's not an issue.)
        os.utime(os.path.join(rosetta_dir,"source",f),None)

        subtest_results = json_results['tests'].setdefault(f, {})
        state = subtest_results.get(_StateKey_, _S_passed_)

        log_line = "{}:{} {}".format(f,l,m)
        if m in historical.get(f,[]):
            historical_log.append(log_line)
        else:
            state = _S_failed_
            problem_log.append(log_line)

        subtest_results[ _StateKey_ ] = state
        subtest_results[ _LogKey_ ] = subtest_results.get( _LogKey_, '' ) + log_line + '\n'
        if issues.get((f,l,m), []): # Has associated html files
            subtest_results.setdefault( _HtmlKey_, [] ).extend( [ p+"#EndPath" for p in issues.get((f,l,m), []) ] )

        if state == _S_failed_ and f not in json_results['summary']['failed_tests']:
            json_results['summary']['failed_tests'].append(f)
    issues_log = ""
    if len(problem_log):
        issues_log += '\nProblems found:\n\n'
        issues_log += '\n'.join(problem_log) + '\n'
    if len(historical_log):
        issues_log += '\nRemaining historical issues:\n\n'
        issues_log += '\n'.join(historical_log) + '\n'

    issues_log += "\nNew issues found: " + str(len(problem_log)) + "\n"
    issues_log += "Total issues found: " + str(len(problem_log) + len(historical_log)) + "\n"
    ## The test server actually has a timer on it.
    #issues_log += "\nRun walltime: " + str(int(time.time()-start_time)) + " s.\n"

    if run_res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = 'Error when running clang static analysis...: {}\nClang Analysis Version: {}\n\n'.format(run_command_line, clang_version) + full_log
        return results

    if len(problem_log) > 0:
        results[_StateKey_] = _S_failed_
    else:
        results[_StateKey_] = _S_passed_

    results[_LogKey_]    =  'Running: {}\nClang Analysis Version: {}\n'.format(run_command_line, clang_version) + '\n' + issues_log
    results[_ResultsKey_] = json_results

    if os.path.exists( os.path.join( working_dir, 'index.html' ) ):
        results[_HtmlKey_] = 'index.html'

    return results



def run_beautification_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Check if branch diff against master is beutified
    '''
    jobs = config['cpu_count']

    execute('Updating options, ResidueTypes and version info...', 'cd {}/source && ./update_options.sh && ./update_ResidueType_enum_files.sh && python version.py'.format(rosetta_dir) )

    state, results, output = _S_script_failed_, {}, ''
    res, o = execute('Reading the source code, trying to see a beauty in it...', 'cd {}/source && python ../../tools/python_cc_reader/test_all_files_already_beautiful.py -j {}'.format(rosetta_dir, jobs), return_='tuple')

    if res:
        state, output = _S_failed_, 'Some of the source code looks ugly!!! Script test_all_files_already_beautiful.py output:\n' + o
    else:
        state = _S_passed_

    return {_StateKey_ : state,  _ResultsKey_ : results,  _LogKey_ : output }



def run_beautify_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Find out to which branch current commit belong to and try to beautify it
    '''
    jobs = config['cpu_count']

    state, results, output = _S_script_failed_, {}, ''

    execute('Updating options, ResidueTypes and version info...', 'cd {}/source && ./update_options.sh && ./update_ResidueType_enum_files.sh && python version.py'.format(rosetta_dir) )

    res, _ = execute('Checking if there is local changes in main repository...', 'cd {} && ( git --no-pager diff --no-color --exit-code >/dev/null || git --no-pager diff --no-color --exit-code --cached >/dev/null ) '.format(rosetta_dir), return_='tuple')
    if res:
        state, output = _S_failed_, 'Working directory is not clean! `git status`:\n{}\nThis might be because you trying to beautify pull-request..., please try too schedule `beautify` test for a branch or SHA1...\n'.format( execute('Checking if there is local changes in main...', 'cd {} && git status'.format(rosetta_dir), return_='output') )

    else:
        o = execute('Getting list of branches where HEAD is present...', 'cd {} && git branch --no-color --contains HEAD && git branch --no-color -r --contains HEAD'.format(rosetta_dir), return_='output')
        branches = set()
        for line in o.split('\n'):
            br = line.replace('*', '').split()
            if br:
                br = br[0]
                if br.startswith('origin/'):
                    br = br[len('origin/'):]
                    #if br != 'master': branches.add(br)
                    branches.add(br)
                #if br != 'HEAD': branches.add(br)

        if 'master' in branches: branches = set( ['master'] )

        commit = execute('Getting current commit sha1...', 'cd {} && git rev-parse HEAD'.format(rosetta_dir), return_='output')
        if len(branches) != 1:
            state, output = _S_failed_, 'Could not figure out which branch to beautify, commit {} belong to following branches (expecting exactly one origin/<branch> or origin/master to be present):\n{}\nAborting...\n'.format(commit, o)

        else:
            branch = branches.pop()
            output += 'Beautifying branch: {} at {} \n'.format(branch, commit)

            output += execute('Checking out branch...', 'cd {} && git fetch && git update-ref refs/heads/{branch} origin/{branch} && git reset --hard {branch} && git checkout {branch} && git branch --set-upstream-to=origin/{branch} {branch}'.format(rosetta_dir, branch=branch), return_='output', add_message_and_command_line_to_output=True)

            if branch == 'master': res, o = execute('Beautifying...', 'cd {}/source && python ../../tools/python_cc_reader/beautify_rosetta.py --overwrite -j {}'.format(rosetta_dir, jobs), return_='tuple', add_message_and_command_line_to_output=True)
            else: res, o = execute('Beautifying...', 'cd {}/source && python ../../tools/python_cc_reader/beautify_changed_files_in_branch.py -j {}'.format(rosetta_dir, jobs), return_='tuple', add_message_and_command_line_to_output=True)

            if res:
                state, output = _S_failed_, 'Beautification script failed with output: {}\nAborting...\n'.format(o)

            else:
                output += o

                res, _ = execute('Checking if there is local changes in main repository...', 'cd {} && ( git --no-pager diff --no-color --exit-code >/dev/null && git --no-pager diff --no-color --exit-code --cached >/dev/null ) '.format(rosetta_dir), return_='tuple')

                if res:
                    output += execute('Calculating changes...', "cd {} && git --no-pager diff --no-color".format(rosetta_dir), return_='output')
                    res, o = execute('Committing and pushing changes...', "cd {} && git commit -a -m 'beautifying' --author='{user_name} <{user_email}>' && git fetch && git rebase && git push".format(rosetta_dir, branch, user_name=config['user_name'], user_email=config['user_email']), return_='tuple')

                    output += o
                    if res: state = _S_failed_

                else: output += '\nBeautification script finished: no beautification required for branch {}!\n'.format(branch)

                output = output.replace( '<{}>'.format(config['user_email']), '<rosetta@university.edu>' )

                state =_S_passed_



    return {_StateKey_ : state,  _ResultsKey_ : results,  _LogKey_ : output }



def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if   test == 'serialization':  return run_serialization_test (rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'clang_analysis': return run_clang_analysis_test(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'beautification': return run_beautification_test(rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'beautify':       return run_beautify_test      (rosetta_dir, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Build script does not support TestSuite-like run!')
