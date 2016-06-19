#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   tests/release.py
## @brief  Rosetta and PyRosetta release scripts
## @author Sergey Lyskov

import os, os.path, json, commands, shutil, tarfile

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'

_number_of_rosetta_binary_revisions_to_keep_in_git_ = 1
_number_of_py_rosetta_revisions_to_keep_in_git_ = 1


def rosetta_source_release(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    TR = Tracer(verbose)
    TR('Running Rosetta source release: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    release_name = 'rosetta.source.{}-{}'.format(config['branch'], config['revision'])
    archive = working_dir + '/' + release_name + '.tar.bz2'

    # Creating git repository with source code, only for regular (not 'commits') branches
    #if config['branch'] != 'commits':
    git_repository_name = 'rosetta.source.{}'.format(config['branch'])
    release_path = '{}/rosetta/git/{}/'.format(config['release_dir'], config['branch'])
    git_origin = os.path.abspath(release_path + git_repository_name + '.git')  # bare repositiry
    git_working_dir = working_dir + '/' + git_repository_name
    if not os.path.isdir(release_path): os.makedirs(release_path)
    if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', 'git init --bare {git_origin} && cd {git_origin} && git update-server-info'.format(**vars()) )

    execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

    # Removing all old files but preserve .git dir...
    execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .*'.format(**vars()), return_='tuple')

    execute('Clonning current checkout of rosetta main...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir} main'.format(**vars()))
    execute('Clonning current checkout of rosetta tools...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../tools tools'.format(**vars()))
    execute('Clonning current checkout of rosetta demos...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../demos demos'.format(**vars()))
    execute('Clonning current checkout of rosetta documentation...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../documentation documentation'.format(**vars()))

    # DANGER DANGER DANGER     DEBUG ONLY, REMOVE LINES BELOW BEFORE COMMITING!!!!!
    #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/convert_to_release.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))
    #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/detect_itest_exes.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))

    execute('Convertion sources to release form...', 'cd {working_dir}/{git_repository_name} && ./tools/release/convert_to_release.bash'.format(**vars()))

    # Creating tar.bz2 archive with sources
    with tarfile.open(archive, "w:bz2") as t: t.add(working_dir+'/'+git_repository_name, arcname=release_name)
    release_path = '{}/rosetta/archive/{}/source/'.format(config['release_dir'], config['branch'])  # , platform['os']
    if not os.path.isdir(release_path): os.makedirs(release_path)

    execute('Moving back upstream .git dir and commiting new release...', 'cd {working_dir}/{git_repository_name} && mv ../.git . && git add * && git ci -a -m "{release_name}"'.format(**vars()))

    execute('Building debug build...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} -j{jobs}'.format(**vars()))  # ignoring extras={extras} because we only test unit test on standard build (not static or MPI etc)
    execute('Building unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} cat=test -j{jobs}'.format(**vars()))  # ignoring extras={extras}
    execute('Building release...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}'.format(**vars()))
    execute('Running unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./test/run.py --compiler={compiler} -j{jobs} --mute all'.format(**vars()))  # ignoring --extras={extras}

    # We moving archive and pushing new revision to upstream only *after* all test runs passed
    shutil.move(archive, release_path+release_name+'.tar.bz2')
    execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --prune=now && git remote prune origin && git push -f'.format(**vars()))

    execute('Pruning origin...', 'cd {git_origin} && git gc --prune=now'.format(**vars()))


    results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : '' }
    json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results


def rosetta_source_and_binary_release(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    TR = Tracer(verbose)
    TR('Running Rosetta source release: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    release_name = 'rosetta.binary.{}.{}-{}'.format(platform['os'], config['branch'], config['revision'])
    archive = working_dir + '/' + release_name + '.tar.bz2'

    # Creating git repository with source code, only for regular (not 'commits') branches
    #if config['branch'] != 'commits':
    git_repository_name = 'rosetta.binary.{}.{}'.format(platform['os'], config['branch'])
    release_path = '{}/rosetta/git/{}/'.format(config['release_dir'], config['branch'])
    git_origin = os.path.abspath(release_path + git_repository_name + '.git')  # bare repositiry
    git_working_dir = working_dir + '/' + git_repository_name
    if not os.path.isdir(release_path): os.makedirs(release_path)
    if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', 'git init --bare {git_origin} && cd {git_origin} && git update-server-info'.format(**vars()) )

    execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

    # Removing all old files but preserve .git dir...
    execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .*'.format(**vars()), return_='tuple')

    execute('Clonning current checkout of rosetta main...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir} main'.format(**vars()))
    execute('Clonning current checkout of rosetta tools...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../tools tools'.format(**vars()))
    execute('Clonning current checkout of rosetta demos...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../demos demos'.format(**vars()))
    execute('Clonning current checkout of rosetta documentation...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../documentation documentation'.format(**vars()))

    # DANGER DANGER DANGER     DEBUG ONLY, REMOVE LINE BELOW BEFORE COMMITING!!!!!
    #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/convert_to_release.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))
    #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/detect_itest_exes.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))

    execute('Convertion sources to release form...', 'cd {working_dir}/{git_repository_name} && ./tools/release/convert_to_release.bash'.format(**vars()))

    execute('Building release...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}'.format(**vars()))

    # Creating tar.bz2 archive with sources
    with tarfile.open(archive, "w:bz2") as t: t.add(working_dir+'/'+git_repository_name, arcname=release_name)
    release_path = '{}/rosetta/archive/{}/binary.{}/'.format(config['release_dir'], config['branch'], platform['os'])
    if not os.path.isdir(release_path): os.makedirs(release_path)

    execute('Moving back upstream .git dir and commiting new release...', 'cd {working_dir}/{git_repository_name} && mv ../.git .'.format(**vars()))
    execute('Adding files and commiting new release...', 'cd {working_dir}/{git_repository_name} && git add * && git add -f main/source/bin/* main/source/build/* && git ci -a -m "{release_name}"'.format(**vars()))

    res, oldest_sha = execute('Getting HEAD~N old commit...', 'cd {working_dir}/{git_repository_name} && git rev-parse HEAD~{_number_of_rosetta_binary_revisions_to_keep_in_git_}'.format(_number_of_rosetta_binary_revisions_to_keep_in_git_=_number_of_rosetta_binary_revisions_to_keep_in_git_, **vars()), return_='tuple')
    if not res:  # if there is no histore error would be raised, but that also mean that rebase is not needed...
        git_truncate = 'git checkout --orphan _temp_ {oldest_sha} && git commit -m "Truncating git history" && git rebase --onto _temp_ {oldest_sha} master && git checkout master && git branch -D _temp_'.format(**vars())
        execute('Trimming git history...', 'cd {working_dir}/{git_repository_name} && {git_truncate}'.format(**vars()))

    # Running extra test to make sure our release is good...
    execute('Building debug build...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} -j{jobs}'.format(**vars()))  # ignoring extras={extras} because we only test unit test on standard build (not static or MPI etc)
    execute('Building unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} cat=test -j{jobs}'.format(**vars()))  # ignoring extras={extras}
    execute('Running unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./test/run.py --compiler={compiler} -j{jobs} --mute all'.format(**vars()))  # ignoring --extras={extras}

    # We moving archive and pushing new revision to upstream only *after* all test runs passed
    shutil.move(archive, release_path+release_name+'.tar.bz2')
    execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --prune=now && git remote prune origin && git push -f'.format(**vars()))

    execute('Pruning origin...', 'cd {git_origin} && git gc --prune=now'.format(**vars()))

    results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : '' }
    json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results


def py_rosetta_release(kind, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread
    #kind = dict(monolith='monolith', namespace='namespace')[ platform['options']['py'] ]  # simple validation: build kind should be ether monolith or namespace,
    #kind_option = '--monolith' if kind == 'monolith' else ''
    kind_option = dict(monolith='--monolith', namespace='', monolith_debug='--monolith --debug', namespace_debug='--debug')[kind]

    TR = Tracer(verbose)

    TR('Running PyRosetta release: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])
    command_line = 'cd {rosetta_dir}/source && BuildPyRosetta.sh -u {kind_option} -j{jobs}'.format(rosetta_dir=rosetta_dir, compiler=compiler, jobs=jobs, extras=extras, kind_option=kind_option)

    if debug: res, output = 0, 'release.py: debug is enabled, skippig build phase...\n'
    else:
        res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line), return_='tuple')
        if res:  res, output = execute('Compiling...', 'cd {}/source && {}'.format(rosetta_dir, command_line.format(compiler=compiler, jobs=1, extras=extras)), return_='tuple')

    file(working_dir+'/build-log.txt', 'w').write(output)

    if res:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)

    else:
        buildings_path_output = execute('Getting buindings build path...', command_line + ' --print-build-path', return_='tuple')
        buildings_path = buildings_path_output[1].split()[-1]
        if not (buildings_path  and  os.path.isdir(buildings_path) ): raise BenchmarkError('Could not retrieve valid PyRosetta bindings binary path!\nCommand line:{}\nResult:{}\n'.format(command_line, buildings_path_output))
        TR('Bindings build path is:{}'.format(buildings_path))

        shutil.copy(config['boost_python_library'], buildings_path)

        memory = config['memory'];  jobs = config['cpu_count']
        if platform['os'] != 'windows': jobs = jobs if memory/jobs >= 2.0 else max(1, int(memory/2) )  # PyRosetta tests require at least 2Gb per memory per thread

        #distr_file_list = os.listdir(buildings_path)

        if debug  or kind.endswith('debug'): res, output = 0, 'release.py: debug is enabled, skippig test phase...\n'
        else:
            res, output = execute('Running PyRosetta tests...', 'cd {buildings_path} && python TestBindings.py -j{jobs}'.format(buildings_path=buildings_path, jobs=jobs), return_='tuple')

        # json_file = buildings_path + '/.test_bindings.json'
        # results = json.load( file(json_file) )

        # execute('Deleting PyRosetta tests output...', 'cd {buildings_path} && python TestBindings.py --delete-tests-output'.format(buildings_path=buildings_path), return_='tuple')
        # extra_files = [f for f in os.listdir(buildings_path) if f not in distr_file_list]  # not f.startswith('.test.')  and
        # if extra_files:
        #     results['results']['tests']['TestBindings'] = dict(state='failed', log='TestBindings.py scripts failed to delete files: ' + ' '.join(extra_files))
        #     results[_StateKey_] = 'failed'

        if not res: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
        output = 'Running: {}\n'.format(command_line) + output  # Making sure that exact command line used is stored

        if res:
            json_file = buildings_path + '/.test.output/.test.results.json'
            results = json.load( file(json_file) )
            results[_LogKey_] = output
            #res_code = _S_build_failed_
            #results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)
        else:
            release_name = 'PyRosetta.{kind}.{os}.{branch}-{revision}'.format(kind=kind, os=platform['os'], branch=config['branch'], revision=config['revision'])
            archive = working_dir + '/' + release_name + '.tar.bz2'

            file_list = 'app database demo test toolbox PyMOLPyRosettaServer.py SetPyRosettaEnvironment.sh TestBindings.py libboost_python rosetta.so'.split()  # rosetta dir is spefial, we omit it here  # ignore_list: _build_ .test.output
            # debug: file_list = 'app demo test toolbox PyMOLPyRosettaServer.py SetPyRosettaEnvironment.sh TestBindings.py libboost_python'.split()  #  ignore_list: _build_ .test.output

            # Creating tar.bz2 archive with binaries
            def arch_filter(tar_info):
                if tar_info.name == release_name: return tar_info
                if tar_info.name == release_name+ '/database': tar_info.type = tarfile.DIRTYPE

                for e in file_list:
                    if tar_info.name.startswith(release_name+'/'+e): return tar_info

                if tar_info.name.startswith(release_name+'/rosetta'):
                    if tar_info.name == release_name+'/rosetta/config.json': return tar_info
                    elif tar_info.type == tarfile.DIRTYPE: return tar_info
                    elif tar_info.name.startswith(release_name+'/rosetta/'):  # special filtering for namespace rosetta dir
                        for ending in '.so .py .pyc .pyd .dylib'.split():
                            if tar_info.name.endswith(ending): return tar_info
                        else: return None

                return None
            with tarfile.open(archive, "w:bz2") as t: t.add(buildings_path, arcname=release_name, filter=arch_filter)

            release_path = '{release_dir}/PyRosetta/archive/{branch}/{kind}.{os}/'.format(release_dir=config['release_dir'], branch=config['branch'], kind=kind, os=platform['os'])
            if not os.path.isdir(release_path): os.makedirs(release_path)
            shutil.move(archive, release_path+release_name+'.tar.bz2')


            # Creating git repository with binaries, only for named branches
            if config['branch'] != 'commits' or True:
                git_repository_name = 'PyRosetta.{kind}.{os}.{branch}'.format(kind=kind, os=platform['os'], branch=config['branch'])
                release_path = '{}/PyRosetta/git/{}/'.format(config['release_dir'], config['branch'])
                git_origin = os.path.abspath(release_path + git_repository_name + '.git')  # bare repositiry
                git_working_dir = working_dir + '/' + git_repository_name
                if not os.path.isdir(release_path): os.makedirs(release_path)
                if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', 'git init --bare {git_origin} && cd {git_origin} && git update-server-info'.format(**vars()) )

                execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

                # Removing all old files but preserve .git dir...
                execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .* ; mv ../.git .'.format(**vars()), return_='tuple')

                for f in os.listdir(buildings_path):
                    for c in file_list:
                        if f.startswith(c):
                            src = buildings_path+'/'+f;  dest = working_dir+'/'+git_repository_name+'/'+f
                            if os.path.isfile(src): shutil.copy(src, dest)
                            elif os.path.isdir(src): shutil.copytree(src, dest)
                            execute('Git add {f}...', 'cd {working_dir}/{git_repository_name} && git add {f}'.format(**vars()))

                build_rosetta_namespace_dir = buildings_path + '/rosetta'
                git_rosetta_namespace_dir = working_dir+'/'+git_repository_name+'/rosetta'
                if os.path.isdir(build_rosetta_namespace_dir):
                    shutil.copytree(build_rosetta_namespace_dir, git_rosetta_namespace_dir)
                    for path, _, files in os.walk(git_rosetta_namespace_dir):
                        for f in files:
                            add = False
                            name = path + '/' + f
                            if name == git_rosetta_namespace_dir+'/config.json': add = True
                            else:
                                for ending in '.so .py .pyc .pyd .dylib'.split():
                                    if name.endswith(ending): add = True

                            if add: execute('Git add {name}...', 'cd {working_dir}/{git_repository_name} && git add {name}'.format(**vars()))


                res, git_output = execute('Git commiting changes...', 'cd {working_dir}/{git_repository_name} && git commit -a -m "{release_name}"'.format(**vars()), return_='tuple')
                if res  and 'nothing to commit, working directory clean' not in git_output: raise BenchmarkError('Could not commit changess to: {}!'.format(git_origin))

                res, oldest_sha = execute('Getting HEAD~N old commit...', 'cd {working_dir}/{git_repository_name} && git rev-parse HEAD~{_number_of_py_rosetta_revisions_to_keep_in_git_}'.format(_number_of_py_rosetta_revisions_to_keep_in_git_=_number_of_py_rosetta_revisions_to_keep_in_git_, **vars()), return_='tuple')
                if not res:  # if there is no histore error would be raised, but that also mean that rebase is not needed...
                    git_truncate = 'git checkout --orphan _temp_ {oldest_sha} && git commit -m "Truncating git history" && git rebase --onto _temp_ {oldest_sha} master && git checkout master && git branch -D _temp_'.format(**vars())  #
                    execute('Trimming git history...', 'cd {working_dir}/{git_repository_name} && {git_truncate}'.format(**vars()))

                execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --prune=now && git remote prune origin && git push -f'.format(**vars()))

                execute('Pruning origin...', 'cd {git_origin} && git gc --prune=now'.format(**vars()))


            #r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            #results[_LogKey_] = output

            res_code = _S_passed_
            results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results




def run(test, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''

    if   test =='source': return rosetta_source_release(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='binary': return rosetta_source_and_binary_release(rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.monolith':  return py_rosetta_release('monolith',  rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.namespace': return py_rosetta_release('namespace', rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.monolith_debug':  return py_rosetta_release('monolith_debug',  rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.namespace_debug': return py_rosetta_release('namespace_debug', rosetta_dir, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    else: raise BenchmarkError('Unknow PyRosetta test: {}!'.format(test))
