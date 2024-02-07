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


_api_version_ = '1.1'

_merge_size_failure_threshold_bytes_ = 1024 * 128


def run_serialization_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    jobs = config['cpu_count']

    setup_for_compile(rosetta_dir)

    executable = install_llvm_tool('serialization_validator', source_location='{}/tools/clang_ast_transform/rosetta-refactor-tool'.format(rosetta_dir), platform=platform, config=config)
    #executable_path = executable.rpartition('/')[0]

    tools_dir = os.path.abspath(rosetta_dir+'/tools')

    json_output_dir = working_dir + '/serialization-validator-output'
    os.mkdir( json_output_dir )

    command_line = 'export PYTHONPATH={tools_dir}/python_cc_reader${{PYTHONPATH+:$PYTHONPATH}} && cd {rosetta_dir}/source && python3 ../tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e '\
    '"python3  {rosetta_dir}/tools/clang_ast_transform/run_serialization_validator_on_file.py --executable_path {executable} --json_output_path {json_output_dir} --filename" -n {jobs}'.format(**vars())

    res, output = execute('Running...', command_line, return_='tuple')
    output = 'Running: ' + command_line + '\n' + output + '\n\n'

    json_files = [f for f in os.listdir(json_output_dir) if f.endswith('.json')]
    if json_files:
        state = _S_failed_ if res else _S_passed_
        if res:
            output += "Non-zero return code {}\n".format(res)
        tests = {}
        results = {_TestsKey_:tests}

        for f in json_files:
            with open(json_output_dir+'/'+f) as _fd: jr = json.load(_fd)
            for k in jr:
                key = k[ len('src/') : ].replace('/', '_')
                if _StateKey_ not in jr[k]  or  jr[k][_StateKey_] != _S_passed_:
                    tests[key] = jr[k]
                    state = _S_failed_
                    output += 'Test {} failed!\n'.format(key)
                else:
                    output += 'Test {} passed!\n'.format(key)

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
    # We need to touch all the broken files we know about later on.

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


def run_clang_tidy_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):

    IGNORE_FAILURES = [ 'src/core/scoring/etable/etrie/TrieCountPair1BC3.cc', #Timeout issues on the test server
                        'src/core/scoring/etable/etrie/TrieCountPair1BC4.cc', #Timeout issues on the test server
                        ]

    # list of Clang-tidy warnings that should not be treated as errors
    WARNING_LINES_TO_IGNORE = [
        #'src/ObjexxFCL/FArray2A.hh:567:3: warning: Returning null reference',
    ]

    CLANG_TIDY_OPTIONS="-quiet -header-filter='.*'"
    CLANG_TIDY_TESTS = ["clang-diagnostic-*",
                        "clang-analyzer-*",
                        "bugprone-*",
                        "misc-*",
                        "modernize-use-nullptr",
                        "modernize-use-override",
                        "-bugprone-forward-declaration-namespace", # Skip - This makes spurious errors when we include the .fwd.hh but not the .hh"
                        "-misc-unconventional-assign-operator", # Skip - we use virtual assignment operators too much for this to be useful
                        "-misc-definitions-in-headers", # Skip - We use a few "include only once" headers (e.g. in the option system) with definitions in them
                        "-misc-misplaced-widening-cast", # Skip - not all that useful, catches a bunch of enum conversions
                        "-clang-analyzer-security.FloatLoopCounter", # Skip - There's a number of places in Rosetta where we attempt to evenly sample across a range.
                        "-clang-analyzer-cplusplus.NewDelete", # Skip - There's issues on the test server with shared_ptr in the system headers.
                        "-clang-analyzer-cplusplus.NewDeleteLeaks", # Ibid.
                            # (These are likely also covered by the ClangSA tests, but with better system header exclusion.)
                        "-clang-diagnostic-c++17-extensions", # Skip - do not warn for using C++17 features without -std=c++17 being enabled.
                        ]
    CLANG_TIDY_TESTS = ','.join( CLANG_TIDY_TESTS )

    jobs = config['cpu_count']
    TR = Tracer(verbose)
    full_log = ''
    results = {}

    TR(f'Running test_suite="{__name__}" at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...')

    run_dir = f'{rosetta_dir}/source/cmake/build_clang_tidy/'

    ver_res, tidy_version = execute('Looking up clang-tidy version...', "clang-tidy --version", return_='tuple' )
    tidy_version = tidy_version.strip()
    if ver_res:
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = "Couldn't figure out clang-tidy version...: \n\n" + tidy_version
        return results

    db_command_line = f'cd {run_dir} && ./build_compile_db.sh'
    db_res, db_output = execute('Building compile db...', db_command_line, return_='tuple' )
    compile_db_file = f'{run_dir}/compile_commands.json'
    if db_res or not os.path.isfile(compile_db_file):
        results[_StateKey_] = _S_build_failed_
        results[_LogKey_]   = "Couldn't create compilation database...: {}\n\n".format(db_command_line) + db_output
        return results

    with open(compile_db_file) as f:
        compile_db = json.load(f)

    prior_run_cache = {}
    prior_run_cache_filename = f'{run_dir}/tidy_results_cache.json'
    prior_run_datestamp = None
    if os.path.isfile( prior_run_cache_filename ):
        prior_run_datestamp = os.path.getmtime( prior_run_cache_filename )
        try:
            with open(prior_run_cache_filename) as f:
                prior_run_cache = json.load(f)
        except:
            print("Prior run cache is malformed.")
            prior_run_cache = {}

        if len(prior_run_cache) and prior_run_cache.get("CLANG_TIDY_TESTS","") != CLANG_TIDY_TESTS:
            print("Clearing cache results due to different tests to run.")
            prior_run_cache = {}
        if len(prior_run_cache) and prior_run_cache.get("CLANG_TIDY_VERSION","") != tidy_version:
            print("Clearing cache results due to different clang tidy version.")
            prior_run_cache = {}
    print("Loaded previous run cache with",len(prior_run_cache)-2,"entries")

    run_cache = {}
    jobslist = {}
    for entry in compile_db:
        fullfilename = entry["file"]
        filename = entry["file"][len(f'{rosetta_dir}/source/'):]
        if not filename.endswith('.cc'): continue # Only bother with our cc files
        if filename.endswith('.cxxtest.cc'): continue # We haven't autogenerated any test files
        if filename.startswith('src/apps/pilot') or filename.startswith('src/devel'): continue
        if filename.endswith('utility/libsvm/Svm.cc'): continue # Effectively an external file

        dependency_cache = f'{run_dir}/dependencies/'+filename+'.d'
        # We just to the loop so we have something to break out of
        rerun_needed = True
        while True:
            # If any of these are true, we can skip to the rerun
            if filename not in prior_run_cache:
                print("Running", filename, "due to it not being in the cache.")
                break
            if prior_run_cache[filename]['result'] != 0:
                print("Rerunning", filename, "due to it not having completed successfully last time.")
                break
            if not os.path.isfile( dependency_cache ):
                print("Rerunning", filename, "due to it not having a dependency file.")
                break
            if os.path.getmtime( dependency_cache ) > prior_run_datestamp:
                print("Rerunning", filename, "as it has changed since the last run.")
                break
            rerun_needed = False
            with open( dependency_cache ) as f:
                for line in f:
                    if len(line) == 0 or line[0] != ' ': continue
                    line = line.split()
                    if len(line) == 0: continue
                    depfile = line[0]
                    if os.path.getmtime( depfile ) > prior_run_datestamp:
                        print("Rerunning", filename, "as it uses a file", depfile, "which has changed since the last run.")
                        rerun_needed = True
                        break
            break # Always break out of the "loop"

        if rerun_needed:
            tidy_command = f"clang-tidy -p . {CLANG_TIDY_OPTIONS} -checks='{CLANG_TIDY_TESTS}' {fullfilename}"
            dep_command = f'mkdir -p $(dirname {dependency_cache}) && ' + entry['command'] + f' -MM -MF {dependency_cache}'
            jobslist[filename] = f'cd {run_dir} && ( {dep_command}; {tidy_command} )'
        else:
            run_cache[filename] = prior_run_cache[filename]

    print("Launching clang-tidy run for", len(jobslist), "files.")

    start_time = time.time()
    raw_results = parallel_execute('clang-tidy', jobslist, rosetta_dir, working_dir, jobs, time=15)
    print("Ran clang-tidy in", int(time.time() - start_time), "seconds")

    # Apparently, clang-tidy doesn't always return non-zero for errors
    for key, res in raw_results.items():
        # Passing results gives a single line of output
        if res['result'] == 0 and len(res['output'].strip().split('\n')) > 1:
            for line in WARNING_LINES_TO_IGNORE:
                if line in res['output']: break
            else: res['result'] = 127

        run_cache[key] = res

    run_cache["CLANG_TIDY_TESTS"] = CLANG_TIDY_TESTS
    run_cache["CLANG_TIDY_VERSION"] = tidy_version
    with open(prior_run_cache_filename,'w') as f:
        json.dump(run_cache,f)

    results[_ResultsKey_] = dict(tests={}, summary=dict(total=0, failed=0, failures_ignored=0, failed_tests=[]))
    for jobname, res in run_cache.items():
        results[_ResultsKey_]['summary']['total'] += 1
        if 'result' in res and res['result'] != 0:
            results[_ResultsKey_]['summary']['failed'] += 1
            if jobname in IGNORE_FAILURES:
                results[_ResultsKey_]['summary']['failures_ignored'] += 1
            results[_ResultsKey_]['summary']['failed_tests'].append(jobname)
            results[_ResultsKey_]['tests'][jobname] = { _StateKey_: _S_failed_, _LogKey_: res['output'] }

    nfailed = results[_ResultsKey_]['summary']['failed'] - results[_ResultsKey_]['summary']['failures_ignored']
    if nfailed == 0:
        results[_StateKey_] = _S_passed_
        results[_LogKey_]   = f"Clang tidy with `-checks={CLANG_TIDY_TESTS}` completed successfully.\n\nClang Tidy Version:\n"+tidy_version
    else:
        results[_StateKey_] = _S_failed_
        results[_LogKey_]   = f"Clang tidy with `-checks={CLANG_TIDY_TESTS}` found errors in {nfailed} files.\n\nClang Tidy Version:\n"+tidy_version

    return results


def run_cppcheck_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):

    jobs = config['cpu_count']
    extras   = ','.join(platform['extras'])

    TR = Tracer(verbose)
    TR('Running test cppcheck at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, hpc_driver={hpc_driver}...'.format( **vars() ) )

    command_line = f'cd src && bash ../../tests/benchmark/util/do_cppcheck.sh -j {jobs} -e "{extras}" -w "{working_dir}"'

    if debug: res, full_output = 0, 'build.py: debug is enabled, skipping build phase...\n'
    else: res, full_output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line), return_='tuple')

    codecs.open(working_dir+'/run-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write( 'Running: {}\n{}\n'.format(command_line, full_output) )

    res_code = _S_failed_ if res else _S_passed_

    # Coallate subtest failures.
    subtest_results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]))
    if res_code == _S_failed_:
        # Only check for failed tests if we fail globally.
        subtest_failures = dict()
        for line in full_output.split('\n'):
            if "historical issues" in line:
                break # Ignore historical issues.
            # Error lines like "[core/scoring/hbonds/hbonds.hh:75]: (performance) Function parameter 'current_path' should be passed by reference. [passedByValue]"
            if '[' not in line or ':' not in line:
                continue # Our finds won't work.
            line = line.strip()
            bracket = line.index('[')
            colon = line.index(':')
            if bracket >= colon:
                continue # Odd order.
            test = line[bracket+1:colon]
            subtest_failures.setdefault( test, [] ).append( line )

        for failed_subtest in subtest_failures:
            subtest_results['tests'][failed_subtest] = { _StateKey_: _S_failed_, _LogKey_: '\n'.join(subtest_failures[failed_subtest]) }
            subtest_results['summary']['failed_tests'].append(failed_subtest)
            subtest_results['summary']['total'] += 1
            subtest_results['summary']['failed'] += 1

    if not res and len(full_output) > 64*1024:
        output = '...truncated...\n'+'\n'.join( full_output.split('\n')[-32:] )  # truncating log for passed builds.
    else:
        # Truncation moved to benchmark.py
        output = full_output

    output = 'Running: {}\n'.format(command_line) + output  # Making sure that exact command line used is stored

    r = {_StateKey_ : res_code,  _ResultsKey_ : subtest_results,  _LogKey_ : output }

    short_r = {_StateKey_ : res_code,  _ResultsKey_ : subtest_results } # Without log
    with open(working_dir+'/output.json', 'w') as f: json.dump(short_r, f, sort_keys=True, indent=2)

    return r


def run_beautification_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Check if branch diff against main is beautified
    '''
    jobs = config['cpu_count']

    setup_for_compile(rosetta_dir)

    state, results, output = _S_script_failed_, {}, ''
    res, o = execute('Reading the source code, trying to see a beauty in it...', 'cd {}/source && python3 ../tools/python_cc_reader/test_all_files_already_beautiful.py -j {}'.format(rosetta_dir, jobs), return_='tuple')

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

    setup_for_compile(rosetta_dir)

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
                    branches.add(br)
                #if br != 'HEAD': branches.add(br)

        if 'main' in branches: branches = set( ['main'] )

        commit = execute('Getting current commit sha1...', 'cd {} && git rev-parse HEAD'.format(rosetta_dir), return_='output')
        if len(branches) != 1:
            state, output = _S_failed_, 'Could not figure out which branch to beautify, commit {} belong to following branches (expecting exactly one origin/<branch> or origin/main to be present):\n{}\nAborting...\n'.format(commit, o)

        else:
            branch = branches.pop()
            output += 'Beautifying branch: {} at {} \n'.format(branch, commit)

            output += execute('Checking out branch...', 'cd {} && git fetch && git update-ref refs/heads/{branch} origin/{branch} && git reset --hard {branch} && git checkout {branch} && git submodule update && git branch --set-upstream-to=origin/{branch} {branch}'.format(rosetta_dir, branch=branch), return_='output', add_message_and_command_line_to_output=True)

            if branch == 'main': res, o = execute('Beautifying...', 'cd {}/source && python3 ../tools/python_cc_reader/beautify_rosetta.py --overwrite -j {}'.format(rosetta_dir, jobs), return_='tuple', add_message_and_command_line_to_output=True)
            else: res, o = execute('Beautifying...', 'cd {}/source && python3 ../tools/python_cc_reader/beautify_changed_files_in_branch.py -j {}'.format(rosetta_dir, jobs), return_='tuple', add_message_and_command_line_to_output=True)

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


def run_submodule_regression_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    jobs = config['cpu_count']

    # These are the Rosetta-originated modules which we're likely interested in tracking the most recent submodule main branch
    # This doesn't include those repos (e.g. external code repos) where we're pinning things to a particular revision for a long time.
    # The value is the "main" branch for each.
    submodules_to_check = {
        'demos':'main',
        'documentation':'master',
        'PyRosetta.notebooks':'master',
        'pyrosetta_scripts':'main',
        'rosetta_scripts_scripts':'main',
        'tools':'main',
        'tests/scientific/data':'main',
        'database/additional_protocol_data':'main',
        'database/protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule':'main',
        'source/external/pybind11':'master',
        'source/src/python/PyRosetta/binder':'master',
    }

    retval = {_StateKey_ : _S_script_failed_,  _ResultsKey_ : {},  _LogKey_ : '' }

    if not os.path.exists(f'{rosetta_dir}/.git'):
        retval[_LogKey_] = "Rosetta is not a git repo -- can't test submodules."
        return retval

    # Git sync to make sure the URLs are up-to-date with the current .gitmodules, and then git fetch to make sure we have the remote's current info.
    res, error_msg = execute('Synchronizing submodules...', f'cd {rosetta_dir} && git submodule sync && cd {submodule} && git fetch', return_='tuple')
    if res: retval[_LogKey_] = "Error syncing submodules: " + error_msg; return retval

    res, rosetta_sha1 = execute('Getting current rosetta SHA1...', f'cd {rosetta_dir} && git rev-parse HEAD', return_='tuple')
    if res: retval[_LogKey_] = "Error getting SHA1: " + rosetta_sha1; return retval
    rosetta_sha1 = rosetta_sha1.strip()

    res, merge_base = execute('Getting merge-base with the rosetta main branch...', f'cd {rosetta_dir} && git merge-base origin/main HEAD', return_='tuple')
    if res: retval[_LogKey_] = "Error getting mergebase: " + merge_base; return retval
    merge_base = merge_base.strip()

    if rosetta_sha1 == merge_base:
        # We're doing the check of main - make sure we don't regress w/r/t both sides
        res, parents = execute("Getting current commit's parents...", f'cd {rosetta_dir} && git rev-list --parents -n 1 HEAD', return_='tuple')
        if res: retval[_LogKey_] = "Error getting parents: " + parents; return retval
        to_check = parents.split()[1:] # First entry is the current commit.
    else:
        to_check = [ merge_base ]

    submodule_states = {}
    for submodule, primary_branch in submodules_to_check.items():
        submodule_status = {_StateKey_ : _S_script_failed_, _LogKey_ : '' }
        retval[_ResultsKey_][ submodule ] = submodule_status
        submodule_states[submodule] = "ERROR WHEN RUNNING"

        if not os.path.exists( f'{rosetta_dir}/{submodule}/.git' ):
            # Git directory in submodule missing means we haven't initialized it, and so we're probably okay with versioning.
            submodule_states[submodule] = "Submodule not initialized"
            del retval[_ResultsKey_][submodule]
            continue

        res, current_sha1 = execute("Getting current submodule commit...", f'cd {rosetta_dir} && git ls-tree HEAD {submodule}', return_='tuple')
        if res: submodule_status[_LogKey_] = f"Error getting submodule SHA1 for submodule {submodule} in HEAD: " + current_sha1; continue
        current_sha1 = current_sha1.split()[2] # format `mode commit SHA1 name`

        ################
        # If we regress against any of the parent commits to_check, we have an issue.
        # (Test against the common anscestors of current commit and main branch)
        regression = False
        for check_sha1 in to_check:
            res, prev_sha1 = execute("Getting parent submodule commit ...", f'cd {rosetta_dir} && git ls-tree {check_sha1} {submodule}', return_='tuple')
            if res: submodule_status[_LogKey_] = f"Error getting submodule SHA1 for submodule {submodule} in {check_sha1}: " + prev_sha1; break
            if len(prev_sha1.strip()) == 0: continue # No submodule on this branch
            prev_sha1 = prev_sha1.split()[2] # format `mode commit SHA1 name`

            if current_sha1 == prev_sha1:
                continue # No regression because it's the same

            res, base = execute("Checking for regression...",f'cd {rosetta_dir}/{submodule} && git merge-base {current_sha1} {prev_sha1}', return_='tuple')
            if res: submodule_status[_LogKey_] = f"Error getting relationship between {current_sha1} and {prev_sha1} in submodule {submodule}: " + base; break
            if base.strip() != prev_sha1: # Current SHA1 does not include the previous one
                regression = True
                break

        if len( submodule_status[_LogKey_] ) != 0: continue # Had an error during running.

        if regression:
            submodule_status[_StateKey_] = _S_script_failed_ # We want the magenta
            submodule_status[_LogKey_] = f"Submodule {submodule} is being reset to state that *removes* commits which were previously in Rosetta"
            submodule_states[submodule] = "INTRODUCES A SUBMODULE REGRESSION"
            continue

        #############
        # Test to see if we're include all of the commits in the submodule which the most recent version of rosetta's main branch has
        # (Test against the most recent version of rosetta's main branch)

        res, main_sha1 = execute("Getting main branch's submodule commit...", f'cd {rosetta_dir} && git ls-tree origin/main {submodule}', return_='tuple')
        if res: submodule_status[_LogKey_] = f"Error getting submodule SHA1 for submodule {submodule} in main: " + main_sha1; continue
        #if len(main_sha1.strip()) == 0: continue # No submodule for main??
        main_sha1 = main_sha1.split()[2] # format `mode commit SHA1 name`

        if main_sha1 != current_sha1:
            res, base = execute("Checking for regression...", f'cd {rosetta_dir}/{submodule} && git merge-base {current_sha1} {main_sha1}', return_='tuple')
            if res: submodule_status[_LogKey_] = f"Error getting relationship between {current_sha1} and rosetta's main branch's version ({main_sha1}) in submodule {submodule}: " + base; continue
            base = base.strip()

            # We want to trigger if there's a divergence between this and the most recent main
            # (This is different from before, which was checking for regression regarding the divergence point of main and main.)
            if base != main_sha1 and base != current_sha1: # We have a submodule which isn't up-to-date with main's main
                submodule_states[submodule] = "NEEDS MERGE"
                submodule_status[_StateKey_] = _S_failed_
                submodule_status[_LogKey_] = f"Submodule {submodule} has commits which diverge with the version in Rosetta's main branch's {submodule}. You need to merge the submodule versions."
                continue

        #############
        # Test to see if we're lagging behind the submodule's primary branch
        # (Test against the submodule's main)

        res, submain_sha1 = execute("Getting submodule primary branch...", f'cd {rosetta_dir}/{submodule} && git rev-parse origin/{primary_branch}', return_='tuple')
        if res:
            submodule_states[submodule] = "NO PRIMARY BRANCH TO TEST!"
            submodule_status[_StateKey_] = _S_failed_
            submodule_status[_LogKey_] = "Error message: " + submain_sha1
            continue

        submain_sha1 = submain_sha1.strip()

        if current_sha1 != submain_sha1:
            res, base = execute("Checking for regression...", f'cd {rosetta_dir}/{submodule} && git merge-base {current_sha1} {submain_sha1}', return_='tuple')
            if res: submodule_status[_LogKey_] = f"Error getting relationship between {current_sha1} and submodule's primary branch head ({submain_sha1}) in submodule {submodule}: " + base; continue
            base = base.strip()

            if base != submain_sha1:
                submodule_states[submodule] = "is not up-to-date with submodule primary branch " + primary_branch
                # This is not an error, just an info message
            else:
                submodule_states[submodule] = 'okay'
        else:
            submodule_states[submodule] = 'okay'

        submodule_status[_StateKey_] = _S_passed_
        submodule_status[_LogKey_] = ''

    # Done all the submodules, now do the summary.
    summary = "\n\n"
    for submodule in submodule_states:
        summary += f"{submodule} -- " + submodule_states[submodule] + '\n'
    retval[_LogKey_] = summary
    #retval[_ResultsKey_] # created in-place
    retval[_StateKey_] = _S_passed_
    for state in retval[_ResultsKey_].values():
        if state[_StateKey_] in [ _S_script_failed_, _S_failed_ ]:
            retval[_StateKey_] = _S_failed_
            break

    return retval



# from Rocco tools/coding_util/find_size_of_branch.py
def pprint_size(size):
    if size < 1024:
        return str(size)+"B"
    size //= 1024
    if size < 1024:
        return str(size)+"KiB"
    size //= 1024
    if size < 1024:
        return str(size)+"MiB"
    size //= 1024
    return str(size)+"GiB"


def run_merge_size_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    merge_size_json_path = f'{working_dir}/find_size_of_branch.json'

    add_on = \
        ( f' --branch {config["merge_head"]}' if config["merge_head"] else '' ) + \
        ( f' --ref {config["merge_base"]}' if config["merge_base"] else '' )

    output = execute('Running...', f'{rosetta_dir}/tools/coding_util/find_size_of_branch.py --size 0 --verbose --diff --json {merge_size_json_path}{add_on}', return_='output')

    with open(f'{working_dir}/find_size_of_branch.output', 'w') as f: f.write(output)

    with open(merge_size_json_path) as f: merge_size = json.load(f)

    total = merge_size['compressed_diff_size']

    if total < _merge_size_failure_threshold_bytes_:
        state = _S_passed_
        prefix = f'Reported merge size: {pprint_size(total)}. This value is below {pprint_size(_merge_size_failure_threshold_bytes_)}, - test passed!'

    elif total < _merge_size_failure_threshold_bytes_ * 16:
        state = _S_failed_
        prefix = f'Reported merge size: {pprint_size(total)}. This value is above threshold: {pprint_size(_merge_size_failure_threshold_bytes_)}, - test failed!'

    else:
        state = _S_script_failed_
        prefix = f'Reported merge size: {pprint_size(total)}. This value is WAY ABOVE threshold: {pprint_size(_merge_size_failure_threshold_bytes_)}, - marking test as script-failed!'

    return {_StateKey_ : state,  _ResultsKey_ : {},  _LogKey_ : prefix + '\n\n' + output }



def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if   test == 'serialization':        return run_serialization_test (repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'clang_analysis':       return run_clang_analysis_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'clang_tidy':           return run_clang_tidy_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'cppcheck':             return run_cppcheck_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'beautification':       return run_beautification_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'beautify':             return run_beautify_test      (repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'submodule_regression': return run_submodule_regression_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test == 'merge_size':           return run_merge_size_test(repository_root, working_dir, platform, config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Build script does not support TestSuite-like run!')
