#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   scripts.py
## @brief  Test suite for Rosetta/PyRosetta scripts
## @author Sergey Lyskov


import os, imp, json, shutil, distutils.dir_util
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'



def run_rosetta_scripts_test(name, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']

    TR = Tracer(verbose)

    jobs = config['cpu_count']

    # Building Rosetta binaries
    res, output, build_command_line = build_rosetta(rosetta_dir, platform, jobs, debug=debug)
    if res: return { _StateKey_ : _S_build_failed_,  _ResultsKey_ : {},
                     _LogKey_ : 'Building rosetta failed!\n{}\n{}\n'.format(build_command_line, output) }

    else:
        script = dict(parse='parse_all_scripts.py', verify='verify_all_scripts_accounted_for.py', validate='validate_all_scripts.py')[name]

        if platform['os'] in 'linux ubuntu': os = 'linux'
        elif platform['os'] == 'mac':        os = 'macos'
        else:                                os = platform['os']

        extras = ','.join(platform['extras']) if platform['extras'] else 'default'

        test_output_json = working_dir + '/' + name + '.json'

        command_line_args = '--rosetta {rosetta_dir} --compiler {compiler} --os {os} --extras {extras} --output-file {test_output_json} --keep-intermediate-files --working-dir {working_dir}'.format(compiler=platform['compiler'], **vars())

        res, output = execute('Running {}...', 'cd  {rosetta_dir}/rosetta_scripts_scripts/tests && python {rosetta_dir}/rosetta_scripts_scripts/tests/{script} {command_line_args}'.format(**vars()), return_='tuple', add_message_and_command_line_to_output=True)

        if res:
            results = {_StateKey_ : _S_script_failed_,  _ResultsKey_ : {},  _LogKey_ : output }

        else:
            with open(test_output_json) as f: json_results = json.load(f)

            results = {_StateKey_ : json_results[_StateKey_],  _ResultsKey_ : { _TestsKey_: json_results['tests']},  _LogKey_ : output + '\n\nTest script log:\n' +json_results['log']}

        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)
        return results




def run_pyrosetta_test(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']

    TR = Tracer(verbose)

    TR('Running scripts.pyrosetta: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    python_environment = get_path_to_python_executable(platform, config)
    ve = setup_python_virtual_environment(working_dir+'/ve', python_environment, 'pytest pytest-json')

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode='MinSizeRel', debug=debug)
    #build_command_line = result.command_line
    #pyrosetta_path = result.pyrosetta_path

    for f in os.listdir(result.pyrosetta_path + '/source'):
        if os.path.islink(result.pyrosetta_path + '/source/' + f): os.remove(result.pyrosetta_path + '/source/' + f)
    distutils.dir_util.copy_tree(result.pyrosetta_path + '/source', working_dir + '/source', update=False)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='replace').write(result.output)

    if result.exitcode:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : result.output }
        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    else:
        output = result.output if result.exitcode else '...\n'+'\n'.join( result.output.split('\n')[-32:] )  # truncating log for passed builds.

        TR('Creating PyRosetta package...')

        # python_version = execute('Getting Python version...', '{python} --version'.format(python=result.python), return_='output').split()[1][:3].replace('.', '')
        # release_name = 'PyRosetta4.{kind}.python{python_version}.{platform}'.format(kind=kind, platform='.'.join([platform['os']]+platform['extras']), python_version=python_version)
        # package_dir = working_dir + '/' + release_name
        package_dir = working_dir + '/package'

        execute('Creating PyRosetta4 distribution package...', '{result.command_line} -s -d --create-package {package_dir}'.format(**vars()))

        TR('Installing PyRosetta into local virtual environment...')
        execute('Installing PyRosetta into local virtual environment...', 'cd {package_dir}/setup && {ve.bin}/python setup.py install'.format(**vars()))
        #install_outout = execute('Installing PyRosetta into local virtual environment...', 'cd {package_dir}/setup && {ve.bin}/python setup.py install'.format(**vars()), return_='output', add_message_and_command_line_to_output=True)
        #output += '\n' + install_outout

        report = working_dir + '/.report.json'

        res, test_output = execute('Running py.test for pyrosetta_scripts repository...', 'cd {rosetta_dir}/pyrosetta_scripts && {ve.bin}/py.test --continue-on-collection-errors -vv --jsonapi --json={report}'.format(**vars()), return_='tuple', add_message_and_command_line_to_output=True)
        codecs.open(working_dir+'/test-log.txt', 'w', encoding='utf-8', errors='replace').write(test_output)

        output += '\n' + test_output

        if os.path.isfile(report):
            json_results = json.load( open(report) )
            with open(report, 'w') as f: json.dump(json_results, f, sort_keys=True, indent=2)

            sub_tests_results = {}
            for a in json_results['included']:
                test = a['attributes']
                name = test['name'].replace('.py::test_', '.')
                if test['outcome'] == 'passed': sub_tests_results[name] = {_StateKey_:_S_passed_, _LogKey_ : '' }
                else: sub_tests_results[name] = {_StateKey_ : _S_failed_,  _LogKey_ : test['call']['longrepr'] }

            res_code = _S_failed_ if res else _S_passed_

        else: res_code, sub_tests_results = _S_script_failed_, {}

        results = {_StateKey_ : res_code,  _ResultsKey_ : { _TestsKey_ : sub_tests_results },  _LogKey_ : output }
        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

        # remove package and virtual-environament dir so they do not use space in Benchmark database
        if not debug:
            shutil.rmtree(package_dir)
            shutil.rmtree(ve.root)

    return results






def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if   test =='rosetta.parse':    return run_rosetta_scripts_test('parse',    rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if   test =='rosetta.verify':   return run_rosetta_scripts_test('verify',   rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    if   test =='rosetta.validate': return run_rosetta_scripts_test('validate', rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='pyrosetta': return run_pyrosetta_test(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    else: raise BenchmarkError('Unknow scripts test: {}!'.format(test))
