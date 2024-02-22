#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   tests/release.py
## @brief  Rosetta and PyRosetta release scripts
## @author Sergey Lyskov

import os, os.path, json, shutil, tarfile, datetime, re as re_module
import codecs

try:
    from setuptools.distutils import dir_util as dir_util_module
except ModuleNotFoundError:
    from distutils import dir_util as dir_util_module


import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.1'

_number_of_rosetta_binary_revisions_to_keep_in_git_ = 1
_number_of_py_rosetta_revisions_to_keep_in_git_ = 1
_number_of_archive_files_to_keep_ = 8
_latest_html_ = 'latest.html'

download_template = '''\
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head>
<title>{distr} {link} Download</title>
<meta http-equiv="REFRESH" content="1; url={link}"></head>
<body></body>
</html>'''

def get_platform_release_name(platform):
    addon = dict(linux='.CentOS', ubuntu='.Ubuntu', mac='', m1='', aarch64='.aarch64.Ubuntu')
    return '.'.join([platform['os']]+platform['extras']) + addon[ platform['os'] ]


def release(name, package_name, package_dir, working_dir, platform, config, release_as_git_repository=True, file=None, use_rosetta_versioning=True):
    ''' Create a release packge: tar.bz2 + git repository
        name - must be a name of what is released without any suffices: rosetta, PyRosetta etc
        package_name - base name for archive (without tar.bz2) that should include name, os, revision, branch + other relevant platform info
        package_dir - location of prepared package
        use_rosetta_versioning - if true resulted tar.bz2-archive will be renamed using Rosetta versioning naming
    '''

    TR = Tracer(True)

    branch = config['branch']
    release_root = config['release_root']

    package_versioning_name = '{package_name}.{branch}-{revision}'.format(package_name=package_name, branch=config['branch'], revision=config['revision'])

    if package_dir:
        TR('Creating tar.bz2 for {name} as {package_versioning_name}...'.format( **vars() ) )
        archive = working_dir + '/' + package_versioning_name + '.tar.bz2'
        with tarfile.open(archive, "w:bz2") as t: t.add(package_dir, arcname=package_versioning_name)  # , filter=arch_filter

    else:
        #assert file.endswith('.tar.bz2')
        archive = file


    if use_rosetta_versioning:
        assert archive.endswith('.tar.bz2')


    release_path = f'{release_root}/{name}/archive/{branch}/{package_name}'
    if not os.path.isdir(release_path): os.makedirs(release_path)

    with FileLock( f'{release_path}/.release.lock' ):
        if use_rosetta_versioning: shutil.move(archive, release_path + '/' + package_versioning_name + '.tar.bz2')
        else: shutil.move(archive, release_path + '/' + os.path.basename(archive) )

        # removing old archives and adjusting _latest_html_
        files = [f for f in os.listdir(release_path) if f != _latest_html_  and  f[0] != '.' ]
        files.sort(key=lambda f: os.path.getmtime(release_path+'/'+f))
        for f in files[:-_number_of_archive_files_to_keep_]: os.remove(release_path+'/'+f)
        if files:
            package_file = files[-1]

            with open(release_path+'/'+_latest_html_, 'w') as h: h.write(download_template.format(distr=name, link=package_file))

            htaccess_file_path = f'{release_path}/.htaccess'

            if os.path.isfile(htaccess_file_path):
                with open(htaccess_file_path) as f: htaccess = f.read()
            else:
                htaccess = ''

            redirect_start = 'RedirectMatch 302 (.*).latest$'
            redirect_line = redirect_start + ' $1' + package_file
            htaccess = re_module.sub(re_module.escape(redirect_start) + '(.*)', redirect_line, htaccess)
            print(f'htaccess: {htaccess}')

            if redirect_line not in htaccess: htaccess += '\n' + redirect_line + '\n'

            with open(htaccess_file_path, 'w') as f: f.write(htaccess)


    # Creating git repository
    if release_as_git_repository:
        TR('Creating git repository for {name} as {package_name}...'.format(**vars()) )

        git_repository_name = '{package_name}.{branch}'.format(**vars())
        git_release_path = f'{config["release_root"]}/{name}/git/{branch}/'
        git_origin = os.path.abspath(git_release_path + git_repository_name + '.git')  # bare repositiry

        git_working_dir = working_dir + '/' + git_repository_name
        if not os.path.isdir(git_release_path): os.makedirs(git_release_path)

        with FileLock(f'{git_release_path}/.{git_repository_name}.release.lock'):

            if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', '( git init --initial-branch master --bare {git_origin} || git init --bare {git_origin} ) && cd {git_origin} && git update-server-info'.format(**vars()) )

            execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

            # Removing all old files but preserve .git dir...
            execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .* ; mv ../.git .'.format(**vars()), return_='tuple')

            for f in os.listdir(package_dir):
                # for c in file_list:
                #     if f.startswith(c):
                src = package_dir+'/'+f;  dest = working_dir+'/'+git_repository_name+'/'+f
                if os.path.isfile(src): shutil.copy(src, dest)
                elif os.path.isdir(src): shutil.copytree(src, dest)
                execute(f'Git add {f}...', f'cd {working_dir}/{git_repository_name} && git add {f}' )

            res, git_output = execute('Git commiting changes...', 'cd {working_dir}/{git_repository_name} && git commit -a -m "{package_name}"'.format(**vars()), return_='tuple')
            if res  and 'nothing to commit, working directory clean' not in git_output: raise BenchmarkError('Could not commit changess to: {}!'.format(git_origin))

            res, oldest_sha = execute('Getting HEAD~N old commit...', 'cd {working_dir}/{git_repository_name} && git rev-parse HEAD~{}'.format(_number_of_py_rosetta_revisions_to_keep_in_git_, **vars()), return_='tuple')
            oldest_sha = oldest_sha.split()[0]  # removing \n at the end of output

            if not res:  # if there is no histore error would be raised, but that also mean that rebase is not needed...
                git_truncate = 'git checkout --orphan _temp_ {oldest_sha} && git commit -m "Truncating git history" && git rebase --onto _temp_ {oldest_sha} master && git checkout master && git branch -D _temp_'.format(**vars())  #
                execute('Trimming git history...', 'cd {working_dir}/{git_repository_name} && {git_truncate}'.format(**vars()))

            #execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --force --prune=now && git remote prune origin && git push -f'.format(**vars()))
            execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git remote prune origin && git push -f'.format(**vars()))

            execute('Pruning origin...', 'cd {git_origin} && git gc --force --prune=now'.format(**vars()))

        if os.path.isdir(git_working_dir): shutil.rmtree(git_working_dir)  # removing git dir to keep size of database small


def convert_to_release(rosetta_dir, working_dir, config, git_repository_name, release_name, tracer):
    ''' Convert Rosetta repostiroty into release mode. This include removing all devel files, checking out submodules and so on...
    '''
    info = generate_version_information(rosetta_dir,
                                        branch=config['branch'],
                                        public_release=True,
                                        revision=config['revision'],
                                        date=datetime.datetime.now(),
                                        package=release_name,
                                        url='https://www.rosettacommons.org',
                                        file_name='{working_dir}/{git_repository_name}/main/.release.json'.format(**vars()))  # we placing this into rosetta/main/ instead of rosetta/ so Rosetta developers could not accidently trigger this unnoticed

    execute('Convertion sources to release form...', 'cd {working_dir}/{git_repository_name} && ./main/tools/release/convert_to_release.bash'.format(**vars()))

    ## These have already been cloned via the convert_to_release.bash script
    #execute('Clonning Binder...', 'cd {working_dir}/{git_repository_name}/main/source/src/python/PyRosetta && git clone https://github.com/RosettaCommons/binder.git && cd binder && git checkout {} && rm -rf .git'.format(info['source']['binder'], **vars()))
    #execute('Clonning Pybind11...', 'cd {working_dir}/{git_repository_name}/main/source/external/ && git clone https://github.com/RosettaCommons/pybind11.git && cd pybind11 && git checkout {} && rm -rf .git'.format(info['source']['pybind11'], **vars()))



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
    release_path = '{}/rosetta/git/{}/'.format(config['release_root'], config['branch'])
    git_origin = os.path.abspath(release_path + git_repository_name + '.git')  # bare repositiry
    git_working_dir = working_dir + '/' + git_repository_name
    if not os.path.isdir(release_path): os.makedirs(release_path)
    if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', 'git init --bare {git_origin} && cd {git_origin} && git update-server-info'.format(**vars()) )

    execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

    # Removing all old files but preserve .git dir...
    execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .*'.format(**vars()), return_='tuple')

    execute('Clonning current checkout of rosetta main...', f'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir} main')
    convert_submodule_urls_from_ssh_to_https(f'{working_dir}/{git_repository_name}/main')
    execute('Initing and updating all submodule...', f'cd {working_dir}/{git_repository_name}/main && git submodule update --init --recursive')

    # execute('Clonning current checkout of rosetta tools...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../tools tools'.format(**vars()))
    # execute('Clonning current checkout of rosetta demos...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../demos demos'.format(**vars()))
    # execute('Clonning current checkout of rosetta documentation...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../documentation documentation'.format(**vars()))

    # DANGER DANGER DANGER     DEBUG ONLY, REMOVE LINES BELOW BEFORE COMMITING!!!!!
    # execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/convert_to_release.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))
    # execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/detect_itest_exes.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))

    #execute('Convertion sources to release form...', 'cd {working_dir}/{git_repository_name} && ./tools/release/convert_to_release.bash'.format(**vars()))
    convert_to_release(rosetta_dir, working_dir, config, git_repository_name, release_name, TR)

    # Creating tar.bz2 archive with sources
    with tarfile.open(archive, "w:bz2") as t: t.add(working_dir+'/'+git_repository_name, arcname=release_name)
    release_path = '{}/rosetta/archive/{}/source/'.format(config['release_root'], config['branch'])  # , platform['os']
    if not os.path.isdir(release_path): os.makedirs(release_path)

    execute('Moving back upstream .git dir and commiting new release...', 'cd {working_dir}/{git_repository_name} && mv ../.git . && git add --force *'.format(**vars()))
    #execute('Adding Binder submodule...', 'cd {working_dir}/{git_repository_name} && git submodule add https://github.com/RosettaCommons/binder.git main/source/src/python/PyRosetta/binder && git submodule update --init --recursive'.format(**vars()))
    #execute('Setting Binder submodule SHA1...', 'cd {working_dir}/{git_repository_name}/main/source/src/python/PyRosetta/binder && git checkout {binder_sha1}'.format(**vars()))
    execute('Commiting new release...', 'cd {working_dir}/{git_repository_name} && git commit -a -m "{release_name}"'.format(**vars()))

    if not debug:
        execute('Building debug build...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} -j{jobs}'.format(**vars()))  # ignoring extras={extras} because we only test unit test on standard build (not static or MPI etc)
        execute('Building unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} cat=test -j{jobs}'.format(**vars()))  # ignoring extras={extras}
        execute('Building release...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}'.format(**vars()))
        execute('Running unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./test/run.py --compiler={compiler} -j{jobs} --mute all'.format(**vars()))  # ignoring --extras={extras}

    # We moving archive and pushing new revision to upstream only *after* all test runs passed
    shutil.move(archive, release_path+release_name+'.tar.bz2')

    # removing old archives and adjusting _latest_html_
    files = [f for f in os.listdir(release_path) if f != _latest_html_]
    files.sort(key=lambda f: os.path.getmtime(release_path+'/'+f))
    for f in files[:-_number_of_archive_files_to_keep_]: os.remove(release_path+'/'+f)
    if files:
        with open(release_path+'/'+_latest_html_, 'w') as h: h.write(download_template.format(distr='rosetta.source', link=files[-1]))

    execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --force --prune=now && git remote prune origin && git push -f'.format(**vars()))

    execute('Pruning origin...', 'cd {git_origin} && git gc --force --prune=now'.format(**vars()))

    results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : '' }
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results



def rosetta_source_and_binary_release(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False, release_as_git_repository=False):
    memory = config['memory'];  jobs = config['cpu_count']
    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    TR = Tracer(verbose)
    TR('Running Rosetta source release: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    release_name = 'rosetta.binary.{}.{}-{}'.format(platform['os'], config['branch'], config['revision'])
    archive = working_dir + '/' + release_name + '.tar.bz2'


    if release_as_git_repository:
        # Creating git repository with source code, only for regular (not 'commits') branches
        #if config['branch'] != 'commits':
        git_repository_name = 'rosetta.binary.{}.{}'.format(platform['os'], config['branch'])
        release_path = '{}/rosetta/git/{}/'.format(config['release_root'], config['branch'])
        git_origin = os.path.abspath(release_path + git_repository_name + '.git')  # bare repositiry
        git_working_dir = working_dir + '/' + git_repository_name
        if not os.path.isdir(release_path): os.makedirs(release_path)
        if not os.path.isdir(git_origin): execute('Origin git repository is not present, initializing...', 'git init --bare {git_origin} && cd {git_origin} && git update-server-info'.format(**vars()) )

        execute('Clonning origin...', 'cd {working_dir} && git clone {git_origin}'.format(**vars()))

        # Removing all old files but preserve .git dir...
        execute('Removing previous files...', 'cd {working_dir}/{git_repository_name} && mv .git .. && rm -r * .*'.format(**vars()), return_='tuple')

        execute('Clonning current checkout of rosetta main...', f'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir} main')
        convert_submodule_urls_from_ssh_to_https(f'{working_dir}/{git_repository_name}/main')
        execute('Initing and updating all submodule...', f'cd {working_dir}/{git_repository_name}/main && git submodule update --init --recursive')

        # execute('Clonning current checkout of rosetta tools...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../tools tools'.format(**vars()))
        # execute('Clonning current checkout of rosetta demos...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../demos demos'.format(**vars()))
        # execute('Clonning current checkout of rosetta documentation...', 'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir}/../documentation documentation'.format(**vars()))

        # DANGER DANGER DANGER     DEBUG ONLY, REMOVE LINE BELOW BEFORE COMMITING!!!!!
        #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/convert_to_release.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))
        #execute('Copying convert_to_release script...', 'cp {rosetta_dir}/../tools/release/detect_itest_exes.bash {working_dir}/{git_repository_name}/tools/release'.format(**vars()))

        #execute('Convertion sources to release form...', 'cd {working_dir}/{git_repository_name} && ./tools/release/convert_to_release.bash'.format(**vars()))
        convert_to_release(rosetta_dir, working_dir, config, git_repository_name, release_name, TR)

        execute('Building release...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}'.format(**vars()))

        # Creating tar.bz2 archive with sources
        with tarfile.open(archive, "w:bz2") as t: t.add(working_dir+'/'+git_repository_name, arcname=release_name)
        release_path = '{}/rosetta/archive/{}/binary.{}/'.format(config['release_root'], config['branch'], platform['os'])
        if not os.path.isdir(release_path): os.makedirs(release_path)

        execute('Moving back upstream .git dir and commiting new release...', 'cd {working_dir}/{git_repository_name} && mv ../.git .'.format(**vars()))
        execute('Adding files and commiting new release...', f'cd {working_dir}/{git_repository_name} && git add * && git add -f main/source/bin/* main/source/build/* && git commit -a -m "{release_name}"')

        res, oldest_sha = execute('Getting HEAD~N old commit...', 'cd {working_dir}/{git_repository_name} && git rev-parse HEAD~{_number_of_rosetta_binary_revisions_to_keep_in_git_}'.format(_number_of_rosetta_binary_revisions_to_keep_in_git_=_number_of_rosetta_binary_revisions_to_keep_in_git_, **vars()), return_='tuple')
        if not res:  # if there is no histore error would be raised, but that also mean that rebase is not needed...
            oldest_sha = oldest_sha.split()[0]
            git_truncate = 'git checkout --orphan _temp_ {oldest_sha} && git commit -m "Truncating git history" && git rebase --onto _temp_ {oldest_sha} master && git checkout master && git branch -D _temp_'.format(**vars())
            execute('Trimming git history...', 'cd {working_dir}/{git_repository_name} && {git_truncate}'.format(**vars()))

        # Running extra test to make sure our release is good...
        execute('Building debug build...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} -j{jobs}'.format(**vars()))  # ignoring extras={extras} because we only test unit test on standard build (not static or MPI etc)
        execute('Building unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} cat=test -j{jobs}'.format(**vars()))  # ignoring extras={extras}
        execute('Running unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./test/run.py --compiler={compiler} -j{jobs} --mute all'.format(**vars()))  # ignoring --extras={extras}

        # We moving archive and pushing new revision to upstream only *after* all test runs passed
        shutil.move(archive, release_path+release_name+'.tar.bz2')
        # removing old archives and adjusting _latest_html_
        files = [f for f in os.listdir(release_path) if f != _latest_html_]
        files.sort(key=lambda f: os.path.getmtime(release_path+'/'+f))
        for f in files[:-_number_of_archive_files_to_keep_]: os.remove(release_path+'/'+f)
        if files:
            with open(release_path+'/'+_latest_html_, 'w') as h: h.write(download_template.format(distr='rosetta.binary', link=files[-1]))

        execute('Pushing changes...', 'cd {working_dir}/{git_repository_name} && git gc --force --prune=now && git remote prune origin && git push -f'.format(**vars()))

        execute('Pruning origin...', 'cd {git_origin} && git gc --force --prune=now'.format(**vars()))

    else:
        git_repository_name = 'rosetta.binary.{}.{}'.format(platform['os'], config['branch'])
        git_working_dir = working_dir + '/' + git_repository_name
        if not os.path.isdir(git_working_dir): os.makedirs(git_working_dir)

        execute('Clonning current checkout of rosetta main...', f'cd {working_dir}/{git_repository_name} && git clone {rosetta_dir} main')
        convert_submodule_urls_from_ssh_to_https(f'{working_dir}/{git_repository_name}/main')
        execute('Initing and updating all submodule...', f'cd {working_dir}/{git_repository_name}/main && git submodule update --init --recursive')

        convert_to_release(rosetta_dir, working_dir, config, git_repository_name, release_name, TR)

        execute('Building release...', f'cd {working_dir}/{git_repository_name}/main/source && ./scons.py bin cxx={compiler} extras={extras} mode=release -j{jobs}')

        # Creating tar.bz2 archive with sources
        with tarfile.open(archive, "w:bz2") as t: t.add(working_dir+'/'+git_repository_name, arcname=release_name)
        release_path = '{}/rosetta/archive/{}/binary.{}/'.format(config['release_root'], config['branch'], platform['os'])
        if not os.path.isdir(release_path): os.makedirs(release_path)

        # Running extra test to make sure our release is good...
        execute('Building debug build...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} -j{jobs}'.format(**vars()))  # ignoring extras={extras} because we only test unit test on standard build (not static or MPI etc)
        execute('Building unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./scons.py cxx={compiler} cat=test -j{jobs}'.format(**vars()))  # ignoring extras={extras}
        execute('Running unit tests...', 'cd {working_dir}/{git_repository_name}/main/source && ./test/run.py --compiler={compiler} -j{jobs} --mute all'.format(**vars()))  # ignoring --extras={extras}

        shutil.rmtree(git_working_dir)

        # We moving archive and pushing new revision to upstream only *after* all test runs passed
        shutil.move(archive, release_path+release_name+'.tar.bz2')
        # removing old archives and adjusting _latest_html_
        files = [f for f in os.listdir(release_path) if f != _latest_html_]
        files.sort(key=lambda f: os.path.getmtime(release_path+'/'+f))
        for f in files[:-_number_of_archive_files_to_keep_]: os.remove(release_path+'/'+f)
        if files:
            with open(release_path+'/'+_latest_html_, 'w') as h: h.write(download_template.format(distr='rosetta.binary', link=files[-1]))


    results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : '' }
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results


def rosetta_documentation(repository_root, working_dir, platform, config, hpc_driver, verbose, debug):
    memory = config['memory'];  jobs = config['cpu_count']
    compiler = platform['compiler']
    extras   = ','.join(platform['extras'])

    print(f'Running Rosetta documentation release: at working_dir={working_dir!r} with rosetta_dir={repository_root}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...')

    release_name = f'rosetta.documentation.{config["branch"]}-{config["revision"]}'
    archive = working_dir + '/' + release_name + '.tar.bz2'

    html_path = f'{repository_root}/source/html'

    html_working_dir = working_dir + '/' + release_name
    if not os.path.isdir(html_working_dir): os.makedirs(html_working_dir)

    #logs = '...'
    if os.path.isdir(html_path): shutil.rmtree(html_path)
    logs = execute('Running Doxygen...', f'cd {repository_root}/source && ./BuildDocs.sh --public', add_message_and_command_line_to_output=True, return_='output')

    dir_util_module.copy_tree(html_path, html_working_dir, update=False)

    with tarfile.open(archive, "w:bz2") as t: t.add(html_working_dir, arcname=release_name)
    release_path = f'{config["release_root"]}/rosetta/archive/{config["branch"]}/documentation'
    if not os.path.isdir(release_path): os.makedirs(release_path)

    shutil.move(archive, release_path + '/' + release_name+'.tar.bz2')
    shutil.rmtree(html_working_dir)

    results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : logs }
    return results



def py_rosetta4_release(kind, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    TR = Tracer(True)

    TR('Running PyRosetta4 release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    # 'release' debug ----------------------------------
    # output = 'dummy\n'


    # python_version = execute('Getting Python version...', '{python} --version'.format(python=platform['python']), return_='output').split()[1][:3].replace('.', '')

    # release_name = 'PyRosetta4.{kind}.python{python_version}.{os}'.format(kind=kind, os=platform['os'], python_version=python_version)
    # package_dir = working_dir + '/' + release_name

    # #execute('Creating PyRosetta4 distribution package...', '{build_command_line} --create-package {package_dir}'.format(**vars()), return_='tuple')
    # dir_util_module.copy_tree('/home/benchmark/rosetta/binder/main/source/build/PyRosetta/linux/clang/pyhton-2.7/minsizerel/build/pyrosetta',
    #                              package_dir, update=False)

    # release('PyRosetta4', release_name, package_dir, working_dir, platform, config)

    # res_code = _S_passed_
    # results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
    # json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, file(working_dir+'/output.json', 'w'), sort_keys=True, indent=2)
    # return results
    # ----------------------------------

    release_name = 'PyRosetta4.{kind}.python{python_version}.{platform}'.format(kind=kind, platform='.'.join([platform['os']]+platform['extras']), python_version=platform['python'].replace('.', '') )

    version_file = working_dir + '/version.json'
    generate_version_information(rosetta_dir, branch=config['branch'], revision=config['revision'], package=release_name, url='http://www.pyrosetta.org', file_name=version_file)  # date=datetime.datetime.now(), avoid setting date and instead use date from Git commit

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode=kind, skip_compile=debug, version=version_file)
    build_command_line = result.command_line
    pyrosetta_path = result.pyrosetta_path

    for f in os.listdir(pyrosetta_path + '/source'):
        if os.path.islink(pyrosetta_path + '/source/' + f): os.remove(pyrosetta_path + '/source/' + f)
    dir_util_module.copy_tree(pyrosetta_path + '/source', working_dir + '/source', update=False)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(result.output)

    if result.exitcode:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : result.output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    else:

        distr_file_list = os.listdir(pyrosetta_path+'/build')

        suite_timeout, test_timeout = (128, 2048) if kind == 'Debug' or platform['os'].startswith('aarch64') else (32, 512)

        #gui_flag = '--enable-gui' if platform['os'] == 'mac' else ''
        gui_flag, res, output = '', result.exitcode, result.output

        #if debug: res, output = 0, 'Release script was invoked with `--debug` flag, - skipping PyRosetta unit tests run...\n'
        if False  and  kind == 'Debug': res, output = 0, 'Debug build, skipping PyRosetta unit tests run...\n'
        else:
            packages = ' '.join( get_required_pyrosetta_python_packages_for_testing(platform) ).replace('>', '=').replace('<', '=')
            python_virtual_environment = setup_persistent_python_virtual_environment(result.python_environment, packages)

            command_line = f'{python_virtual_environment.activate} && cd {result.pyrosetta_path}/build && {python_virtual_environment.python} {rosetta_dir}/source/test/timelimit.py {suite_timeout} {python_virtual_environment.python} self-test.py {gui_flag} -j{jobs} --timeout {test_timeout}'
            output += '\nRunning PyRosetta tests: ' + command_line + '\n'

            res, o = execute('Running PyRosetta tests...', command_line, return_='tuple')
            output += o

            #res, output = execute('Running PyRosetta tests...', 'cd {pyrosetta_path}/build && {python} self-test.py {gui_flag} -j{jobs} --timeout {timeout}'.format(pyrosetta_path=pyrosetta_path, python=result.python, jobs=jobs, gui_flag=gui_flag, timeout=timeout), return_='tuple')


        json_file = pyrosetta_path + '/build/.test.output/.test.results.json'
        with open(json_file) as f: results = json.load(f)

        execute('Deleting PyRosetta tests output...', 'cd {pyrosetta_path}/build && {python} self-test.py --delete-tests-output'.format(pyrosetta_path=pyrosetta_path, python=result.python), return_='tuple')
        extra_files = [f for f in os.listdir(pyrosetta_path+'/build') if f not in distr_file_list]  # not f.startswith('.test.')  and
        if extra_files:
            results['results']['tests']['self-test'] = dict(state='failed', log='self-test.py scripts failed to delete files: ' + ' '.join(extra_files))
            results[_StateKey_] = 'failed'

        if results[_StateKey_] == _S_passed_: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
        output = 'Running: {}\n'.format(build_command_line) + output  # Making sure that exact command line used is stored

        #r = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        results[_LogKey_] = output

        if results[_StateKey_] == _S_failed_:
            with open(working_dir+'/output.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

            # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            #with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

        else:

            TR('Running PyRosetta4 release test: Build and Unit tests passged! Now creating PyRosetta package...')

            package_dir = working_dir + '/' + release_name

            execute('Creating PyRosetta4 distribution package...', '{build_command_line} -sd --create-package {package_dir}'.format(**vars()))

            release('PyRosetta4', release_name, package_dir, working_dir, platform, config, release_as_git_repository = True if kind in [] else False )

            # releasing PyMOL-RosettaServer scripts
            release('PyMOL-RosettaServer', 'PyMOL-RosettaServer.python2',        package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=f'{package_dir}/PyMOL-RosettaServer.py',                use_rosetta_versioning=False)
            release('PyMOL-RosettaServer', 'PyMOL-RosettaServer.python3',        package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=f'{package_dir}/PyMOL-RosettaServer.python3.py',        use_rosetta_versioning=False)
            release('PyMOL-RosettaServer', 'PyMOL-Rosetta-relay-client.python3', package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=f'{package_dir}/PyMOL-Rosetta-relay-client.python3.py', use_rosetta_versioning=False)
            release('PyMOL-RosettaServer', 'PyMOL-Rosetta-relay-client.python2', package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=f'{package_dir}/PyMOL-Rosetta-relay-client.python2.py', use_rosetta_versioning=False)

            # building and releaseing Wheel archive
            #if (platform['python'][0] == '2' or platform['python'] == '3.5')  and  platform['os'] == 'mac': pass
            if platform['python'][0] == '2': pass
            else:
                whell_environment = setup_persistent_python_virtual_environment(result.python_environment, 'setuptools wheel')

                execute('Creating PyRosetta4 distribution Wheel package...', f'{whell_environment.activate} && cd {package_dir}/setup && python setup.py sdist bdist_wheel')
                wheel_file_name = [ f for f in os.listdir( f'{package_dir}/setup/dist' ) if f.endswith('.whl') ][0]
                release('PyRosetta4', release_name+'.wheel', package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=f'{package_dir}/setup/dist/{wheel_file_name}', use_rosetta_versioning=False)

            if os.path.isdir(package_dir): shutil.rmtree(package_dir)  # removing package to keep size of database small

            results[_StateKey_] = _S_passed_
            results[_LogKey_]   = output

            #with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            with open(working_dir+'/output.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results



def py_rosetta4_documentation(kind, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    #if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    TR = Tracer(True)

    TR('Running PyRosetta4-documentation release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    release_root = config['release_root']

    package_name = 'PyRosetta4.{kind}.python{python_version}.{platform}'.format(kind=kind, platform='.'.join([platform['os']]+platform['extras']), python_version=platform['python'].replace('.', ''))
    package_name = '{package_name}.{branch}-{revision}'.format(package_name=package_name, branch=config['branch'], revision=config['revision'])

    version_file = working_dir + '/version.json'
    generate_version_information(rosetta_dir, branch=config['branch'], revision=config['revision'], package=package_name, url='http://www.pyrosetta.org', file_name=version_file)   # date=datetime.datetime.now(), avoid setting date and instead use date from Git commit

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode=kind, skip_compile=debug, version=version_file)

    packages = ' '.join( get_required_pyrosetta_python_packages_for_testing(platform) ).replace('>', '=').replace('<', '=') + ' sphinx==5.2.3'
    python_virtual_environment = setup_persistent_python_virtual_environment(result.python_environment, packages)

    res = result.exitcode
    output = result.output
    build_command_line = result.command_line
    pyrosetta_path = result.pyrosetta_path

    for f in os.listdir(pyrosetta_path + '/source'):
        if os.path.islink(pyrosetta_path + '/source/' + f): os.remove(pyrosetta_path + '/source/' + f)
    dir_util_module.copy_tree(pyrosetta_path + '/source', working_dir + '/source', update=False)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(output)

    if res:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    else:
        documentation_dir = os.path.abspath(working_dir+'/documentation')

        res, output2 = execute('Generating PyRosetta-4 documentation...', "{python_virtual_environment.activate} && {build_command_line} -s -d --documentation {documentation_dir}".format(**vars()), return_='tuple', add_message_and_command_line_to_output=True)

        if res:
            res_code = _S_build_failed_
            results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output+output2 }
            with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

        else:

            with FileLock( f'{release_root}/PyRosetta4/documentation/.release.lock' ):

                release_path = '{release_root}/PyRosetta4/documentation/PyRosetta-4.documentation.{branch}.{kind}.python{python_version}.{os}'.format(branch=config['branch'], os=platform['os'], python_version=platform['python'].replace('.', ''), **vars())

                if os.path.isdir(release_path): shutil.rmtree(release_path)
                shutil.move(documentation_dir, release_path)

                res_code = _S_passed_
                results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output+output2 }
                with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results


_index_html_template_ = '''\
<html>
<head>
    <title>PyRosetta conda package</title>

  <style>
    fixed {{background-color: #eee; white-space: pre-wrap; font-family: Monaco, 'Liberation Mono', Courier, monospace; font-size:12px; }}
  </style>
</head>
<body>
<p>
    To install this PyRosetta conda package:
    <ul>
        <li>
        please add <fixed>graylab.jhu.edu/download/PyRosetta4/conda/{release_kind}</fixed> into your local <fixed>~/.condarc</fixed> file, like:<br/><br/>
<fixed>channels:
  - https://USERNAME:PASSWORD@graylab.jhu.edu/download/PyRosetta4/conda/{release_kind}
  - defaults
</fixed>
<br/><br/>(ask for user-name and password in RosettaCommons Slack <fixed>#PyRosetta</fixed> channel)
        </li>
        <li> Then run <fixed>conda install pyrosetta={conda_package_version}</fixed> to install <em>this<em> build.
        </li>

    </ul>
</p>

</body></html>
'''


_conda_setup_only_build_sh_template_ = '''\
#Configure!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -euo pipefail
IFS=$'\n\t'

set -x

echo "--- Build"
echo "PWD: `pwd`"
echo "Python: `which python` --> `python --version`"
echo "PREFIX Python: `which ${{PREFIX}}/bin/python` --> `${{PREFIX}}/bin/python --version`"


echo "-------------------------------- Installing PyRosetta Python package..."

pushd {package_dir}/setup

cat ../version.json

# Run initial test to prebuild databases
${{PREFIX}}/bin/python -c 'import pyrosetta; pyrosetta.init(); pyrosetta.get_score_function()(pyrosetta.pose_from_sequence("TEST"))'

${{PREFIX}}/bin/python setup.py install --single-version-externally-managed --record=record.txt > install.log

popd
echo "-------------------------------- Installing PyRosetta Python package... Done."
'''

def native_libc_py_rosetta4_conda_release(kind, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    if 'cxx11thread' not in platform['extras']  or  'serialization' not in platform['extras']: raise BenchmarkError( f'Running native_libc_py_rosetta4_conda_release: on platform with extras={platform["extras"]}, however Conda build on platform without cxx11thread or serialization is not supported!' )

    TR = Tracer(True)

    TR('Running PyRosetta4 conda release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    conda = setup_conda_virtual_environment(working_dir, platform, config, packages='setuptools')

    platform_name = get_platform_release_name(platform)
    release_name = 'PyRosetta4.conda.{platform}.python{python_version}.{kind}'.format(kind=kind, platform=platform_name, python_version=platform['python'].replace('.', '') )

    version_file = working_dir + '/version.json'
    version = generate_version_information(rosetta_dir, branch=config['branch'], revision=config['revision'], package=release_name, url='http://www.pyrosetta.org', file_name=version_file)  # date=datetime.datetime.now(), avoid setting date and instead use date from Git commit

    result = build_pyrosetta(rosetta_dir, platform, jobs, config, mode=kind, conda=conda, skip_compile=debug, version=version_file, options='--no-strip-module --binder-config rosetta.distributed.config')  # --multi-threaded  --serialization
    build_command_line = result.command_line
    pyrosetta_path = result.pyrosetta_path

    for f in os.listdir(pyrosetta_path + '/source'):
        if os.path.islink(pyrosetta_path + '/source/' + f): os.remove(pyrosetta_path + '/source/' + f)
    dir_util_module.copy_tree(pyrosetta_path + '/source', working_dir + '/source', update=False)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(result.output)

    if result.exitcode:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : result.output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    else:

        if debug:
            res, output = 0, 'Benchmark `debug` is enabled, skipping PyRosetta unit test run...\n'
            results = {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : output }

        else:
            distr_file_list = os.listdir(pyrosetta_path+'/build')

            packages = ' '.join( get_required_pyrosetta_python_packages_for_testing(platform) ).replace('>', '=').replace('<', '=')

            local_python = local_python_install(platform, config)
            python_virtual_environment = setup_persistent_python_virtual_environment(local_python, packages)

            #gui_flag = '--enable-gui' if platform['os'] == 'mac' else ''
            gui_flag, res, output = '', result.exitcode, result.output
            if False  and  kind == 'Debug': res, output = 0, 'Debug build, skipping PyRosetta unit tests run...\n'
            else: res, output = execute('Running PyRosetta tests...', f'{python_virtual_environment.activate} && cd {pyrosetta_path}/build && {python_virtual_environment.python} {rosetta_dir}/source/test/timelimit.py 128 {python_virtual_environment.python} self-test.py --timeout 2048 {gui_flag} -j{jobs}', return_='tuple')

            json_file = pyrosetta_path + '/build/.test.output/.test.results.json'
            with open(json_file) as f: results = json.load(f)

            execute('Deleting PyRosetta tests output...', 'cd {pyrosetta_path}/build && {python} self-test.py --delete-tests-output'.format(pyrosetta_path=pyrosetta_path, python=result.python), return_='tuple')
            extra_files = [f for f in os.listdir(pyrosetta_path+'/build') if f not in distr_file_list]  # not f.startswith('.test.')  and
            if extra_files:
                results['results']['tests']['self-test'] = dict(state='failed', log='self-test.py scripts failed to delete files: ' + ' '.join(extra_files))
                results[_StateKey_] = 'failed'

            if results[_StateKey_] == _S_passed_: output = '...\n'+'\n'.join( output.split('\n')[-32:] )  # truncating log for passed builds.
            output = 'Running: {}\n'.format(build_command_line) + output  # Making sure that exact command line used is stored

            results[_LogKey_] = output

        if results[_StateKey_] == _S_failed_:
            # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)
        else:

            TR('Running PyRosetta4 release test: Build and Unit tests passged! Now creating PyRosetta package...')

            package_dir = working_dir + '/' + release_name

            execute( f'Creating PyRosetta4 distribution package...', f'{build_command_line} -sd --create-package {package_dir}' )

            python_version_as_tuple = tuple( map(int, platform.get('python', DEFAULT_PYTHON_VERSION).split('.') ) )

            recipe_dir = working_dir + '/recipe';  os.makedirs(recipe_dir)

            recipe = dict(
                package = dict(
                    name    = 'pyrosetta',
                    version = version['version'],
                ),
                requirements = dict(
                    build = [f'python {platform["python"]}'],
                    host  = [f'python {platform["python"]}', 'setuptools', 'zlib'],
                    #run   = [f'python =={platform["python"]}', "{{ pin_compatible('numpy') }}", 'zlib', 'pandas >=0.18', 'scipy >=1.0', 'traitlets', 'python-blosc'],
                    run   = [f'python {platform["python"]}', 'zlib', ] + get_required_pyrosetta_python_packages_for_release_package(platform, conda=True),
                ),

                about = dict(
                    home    = 'http://www.pyrosetta.org',
                    license = 'Rosetta license',
                    license_file = f'{rosetta_dir}/LICENSE.md',
                    summary = 'Python binding for Rosetta, biomolecular modeling software package',
                    description = 'The Rosetta software suite includes algorithms for computational modeling and analysis of protein structures. It has enabled notable scientific advances in computational biology, including de novo protein design, enzyme design, ligand docking, and structure prediction of biological macromolecules and macromolecular complexes.',
                ),
            )

            # now tested in T900_distributed.py
            # if python_version_as_tuple < (3, 9):
            #     recipe['test'] = dict(
            #         requires = [f'python {platform["python"]}'],
            #         commands = ['python -m unittest pyrosetta.tests.distributed.test_smoke']
            #     )

            with open( recipe_dir + '/meta.yaml', 'w' ) as f: json.dump(recipe, f, sort_keys=True, indent=2)

            with open( recipe_dir + '/build.sh', 'w' ) as f: f.write( _conda_setup_only_build_sh_template_.format(**locals()) )

            # --output              Output the conda package filename which would have been created
            # --output-folder OUTPUT_FOLDER folder to dump output package to. Package are moved here if build or test succeeds. Destination folder must exist prior to using this.

            release_kind = 'release' if config['branch'] == 'release' else 'devel'
            conda_release_path = '{release_dir}/PyRosetta4/conda/{release_kind}'.format(release_dir=config['release_root'], release_kind = release_kind)
            if not os.path.isdir(conda_release_path): os.makedirs(conda_release_path)

            with FileLock( '{conda_release_path}/.{os}.python{python_version}.release.lock'.format(os=platform['os'], python_version=platform['python'].replace('.', ''), **vars()) ):
                working_dir_release_path = f'{working_dir}/conda-release'
                os.makedirs(working_dir_release_path)

                conda_build_command_line = f'{conda.activate_base} && conda build purge && conda build --no-locking --quiet {recipe_dir} --output-folder {working_dir_release_path}' # --channel conda-forge
                conda_package_output = execute('Getting Conda package name...', f'{conda_build_command_line} --output', return_='output', silent=True)

                m = re_module.search(r"pyrosetta-.*\.tar\.bz2", conda_package_output, re_module.MULTILINE)
                conda_package = m.group(0) if m else 'unknown'
                conda_package_dir = re_module.search(r"/([^/]*)/pyrosetta-.*\.tar\.bz2", conda_package_output, re_module.MULTILINE).group(1)

                TR(f'Building Conda package: {conda_package}...')
                res, conda_log = execute('Creating Conda package...', conda_build_command_line, return_='tuple', add_message_and_command_line_to_output=True)

                results[_LogKey_]  += f'Got package name from conda build command line `{conda_build_command_line}` : {conda_package}\n' + conda_log
                with open(working_dir+'/conda-build-log.txt', 'w') as f: f.write( to_unicode(conda_log) )

                if not os.path.isdir(f'{conda_release_path}/{conda_package_dir}'): os.makedirs(f'{conda_release_path}/{conda_package_dir}')
                shutil.move(f'{working_dir_release_path}/{conda_package_dir}/{conda_package}', f'{conda_release_path}/{conda_package_dir}/{conda_package}')

            if res:
                results[_StateKey_] = _S_script_failed_
                results[_LogKey_]  += conda_log
            else:
                with FileLock( f'{conda_release_path}/.release.lock' ):
                    execute('Regenerating Conda package index...', f'{conda.activate_base} && cd {conda_release_path} && conda index .')

                conda_package_version = conda_package.split('/')[-1].partition('-')[2]
                print(f'Determent conda_package_version for package {conda_package!r}: {conda_package_version!r}...')
                with open(f'{working_dir}/index.html', 'w') as f: f.write( _index_html_template_.format(**vars() ) )


            if not debug:
                for d in [conda.root, package_dir, working_dir_release_path]: shutil.rmtree(d)  # removing packages to keep size of Benchmark database small

            # res_code = _S_passed_
            # results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            #with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            with open(working_dir+'/output.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space


            '''
            conda_token = config.get('conda_token', '')
            maybe_upload = f' --channel rosettacommons --user rosettacommons --token {conda_token}' if conda_token and config['branch'] == 'release' else ''
            #maybe_upload = f' --channel rosettacommons --user rosettacommons --token {conda_token} --label devel' if conda_token else ''

            conda_build_command_line = f'{conda.activate_base} && conda build purge && conda build --no-locking --quiet {recipe_dir}  --output-folder {conda_package_dir}' + maybe_upload
            conda_package = execute('Getting Conda package name...', f'{conda_build_command_line} --output', return_='output', silent=True).split()[0]  # removing '\n' at the end

            TR(f'Building Conda package: {conda_package}...')
            res, conda_log = execute('Creating Conda package...', conda_build_command_line, return_='tuple', add_message_and_command_line_to_output=True)
            conda_log = conda_log.replace(conda_token, 'CONDA_TOKEN')

            results[_LogKey_]  += f'Got package name from conda build command line `{conda_build_command_line}` : {conda_package}\n' + conda_log
            with open(working_dir+'/conda-build-log.txt', 'w') as f: f.write( to_unicode(conda_log) )

            if res:
                results[_StateKey_] = _S_script_failed_
                results[_LogKey_]  += conda_log

            else:
                release('PyRosetta4', release_name, None, working_dir, platform, config, release_as_git_repository = False, file = conda_package)

            if not debug:
                for d in [conda.root, package_dir, conda_package_dir]: shutil.rmtree(d)  # removing packages to keep size of Benchmark database small

            # res_code = _S_passed_
            # results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
            with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
            '''

    return results


_conda_build_sh_template_ = '''\
#Configure!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -euo pipefail
IFS=$'\n\t'

set -x

echo "--- Build"
echo "PWD: `pwd`"
echo "GCC: `which gcc`"
echo "Python: `which python` --> `python --version`"
echo "PREFIX Python: `which ${{PREFIX}}/bin/python` --> `${{PREFIX}}/bin/python --version`"


echo "--- Env"
build_args=(
--create-package {package_dir}
--version {working_dir}/version.json
--binder-config rosetta.config
--binder-config rosetta.distributed.config
--serialization
--multi-threaded
--no-strip-module
)
#--no-zmq

if [[ ! -z "${{GCC:-}}" ]]; then
  # Build via gcc/g++ rather than conda cc c++ compiler aliases
  # binder invokation still targets system C++ standard library
  # see linux-anvil for system gcc/g++ installation
  build_args+=(--compiler ${{GCC}})
  export CC=${{GCC}}
  export CXX=${{GXX}}

  # Override flags to just include prefix
  export CFLAGS="-I${{PREFIX}}/include"
  export CXXFLAGS="-I${{PREFIX}}/include"
fi

if [[ ! -z "${{CLANG:-}}" ]]; then
  # override conda-provided clang compiler with system clang
  # still links against conda libc++ shared libraries
  export CLANG=/usr/bin/clang
  export CC=${{CLANG}}
  export CLANGXX=/usr/bin/clang++
  export CXX=${{CLANGXX}}
  build_args+=(--compiler ${{CLANG}})
fi


pushd {rosetta_dir}/source/src/python/PyRosetta

${{PREFIX}}/bin/python build.py ${{build_args[@]}} -j{jobs}

echo "-------------------------------- Running PyRosetta unit tests..."
pushd `${{PREFIX}}/bin/python build.py ${{build_args[@]}} --print-build-root`/build
${{PREFIX}}/bin/python self-test.py -j{jobs}
${{PREFIX}}/bin/python self-test.py --delete-tests-output

echo "-------------------------------- Running PyRosetta unit tests... Done."

popd
popd


echo "-------------------------------- Installing PyRosetta Python package..."

pushd {package_dir}/setup

cat ../version.json

# Run initial test to prebuild databases
${{PREFIX}}/bin/python -c 'import pyrosetta; pyrosetta.init(); pyrosetta.get_score_function()(pyrosetta.pose_from_sequence("TEST"))'

${{PREFIX}}/bin/python setup.py install --single-version-externally-managed --record=record.txt > install.log

popd
echo "-------------------------------- Installing PyRosetta Python package... Done."
'''
def conda_libc_py_rosetta4_conda_release(kind, rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Build PyRosetta package using Conda build tools so it will be linked to Conda provided libc
    '''
    memory = config['memory'];  jobs = config['cpu_count']
    if platform['os'] != 'windows': jobs = jobs if memory/jobs >= PyRosetta_unix_memory_requirement_per_cpu else max(1, int(memory/PyRosetta_unix_memory_requirement_per_cpu) )  # PyRosetta require at least X Gb per memory per thread

    TR = Tracer(True)

    TR('Running PyRosetta4 conda release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={jobs}, memory={memory}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    conda = setup_conda_virtual_environment(working_dir, platform, config, packages='gcc')  # gcc cmake ninja

    release_name = 'PyRosetta4.conda.{kind}.python{python_version}.{platform}'.format(kind=kind, platform='.'.join([platform['os']]+platform['extras']), python_version=platform['python'].replace('.', '') )

    version_file = working_dir + '/version.json'
    version = generate_version_information(rosetta_dir, branch=config['branch'], revision=config['revision'], package=release_name, url='http://www.pyrosetta.org', file_name=version_file)  # date=datetime.datetime.now(), avoid setting date and instead use date from Git commit

    python_version_as_tuple = tuple( map(int, platform.get('python', DEFAULT_PYTHON_VERSION).split('.') ) )

    recipe_dir = working_dir + '/recipe';  os.makedirs(recipe_dir)

    recipe = dict(
        package = dict(
            name    = 'pyrosetta',
            version = version['version'],
        ),
        requirements = dict(
            build = [f'python {platform["python"]}', 'gcc', "{{ compiler('c') }}", "{{ compiler('cxx') }}", ], # 'cmake', 'ninja'
            host  = [f'python {platform["python"]}', 'setuptools', 'numpy', 'zlib'],
            run   = [f'python {platform["python"]}', "{{ pin_compatible('numpy') }}", 'zlib', 'pandas >=0.18', 'scipy >=1.0', 'traitlets', 'python-blosc'],
        ),

        about = dict( home ='http://www.pyrosetta.org' ),
    )

    if python_version_as_tuple < (3, 9): recipe['test'] = dict( commands = ['python -m unittest pyrosetta.tests.distributed.test_smoke'] )

    with open( recipe_dir + '/meta.yaml', 'w' ) as f: json.dump(recipe, f, sort_keys=True, indent=2)

    package_dir = working_dir + '/' + release_name

    with open( recipe_dir + '/build.sh', 'w' ) as f: f.write( _conda_build_sh_template_.format(**locals()) )

    # --output              Output the conda package filename which would have been created
    # --output-folder OUTPUT_FOLDER folder to dump output package to. Package are moved here if build or test succeeds. Destination folder must exist prior to using this.

    conda_package_dir = working_dir + '/conda_package';  os.makedirs(conda_package_dir)

    conda_build_command_line = f'{conda.activate_base} && conda build purge && conda build --no-locking --quiet {recipe_dir} --output-folder {conda_package_dir}'
    conda_package = execute('Getting Conda package name...', f'{conda_build_command_line} --output', return_='output', silent=True).split()[0]  # removing '\n' at the end

    TR(f'Building Conda package: {conda_package}...')
    conda_log = execute('Creating Conda package...', conda_build_command_line, return_='output')
    results[_LogKey_]  += conda_log
    with open(working_dir+'/conda-build-log.txt', 'w') as f: f.write( to_unicode(result.output) )


    '''
    release('PyRosetta4', release_name, package_dir, working_dir, platform, config, release_as_git_repository = True if kind in ['Release', 'MinSizeRel'] else False )

    if os.path.isdir(package_dir): shutil.rmtree(package_dir)  # removing package to keep size of database small

    res_code = _S_passed_
    results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space
    '''

    return results




def ui_release(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    TR = Tracer(True)
    TR('Running ui_release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={config[cpu_count]}, memory={config[memory]}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    platform_suffix = platform_to_pretty_string(platform)
    build_path = '{rosetta_dir}/source/build/ui.{platform_suffix}.static'.format(**vars())
    qt_extras = '-spec linux-clang ' if (platform['compiler'] == 'clang' and platform['os'] == 'linux') else ''

    command_line = 'cd {rosetta_dir}/source/src/ui && python update_ui_project.py && cd ../../build && mkdir -p {build_path} && cd {build_path} && {config[qmake.static]} -r ../qt/qt.pro {qt_extras}&& make -j{config[cpu_count]}'.format(**vars())

    if debug: res, output = 0, 'build.py: debug is enabled, skippig build phase...\n'
    else: res, output = execute('Compiling...', command_line, return_='tuple', add_message_and_command_line_to_output=True)

    codecs.open(working_dir+'/build-log.txt', 'w', encoding='utf-8', errors='backslashreplace').write(output)

    if res:
        res_code = _S_build_failed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)

    else:
        does_not_require_database = ''.split()  # bundle_gui

        apps = 'workbench parametric_design rna_denovo pose_viewer'.split()
        for a in apps:
            release_name = 'ui.{a}.{platform}'.format(a=a, platform='.'.join([platform['os']]) ) #, python_version=platform['python'].replace('.', '') )
            package_dir = working_dir + '/' + release_name
            if not os.path.isdir(package_dir): os.makedirs(package_dir)

            if platform['os'] == 'mac':
                dir_util_module.copy_tree(build_path + '/{a}/{a}.app'.format(**vars()), '{package_dir}/{a}.app'.format(**vars()), update=False)
                if a not in does_not_require_database: dir_util_module.copy_tree(rosetta_dir + '/database', '{package_dir}/{a}.app/Contents/database'.format(**vars()), update=False)

            elif platform['os'] in ['linux', 'ubuntu']:
                #os.makedirs( '{package_dir}/{a}'.format(**vars()) )
                #shutil.copy(build_path + '/{a}/{a}'.format(**vars()), package_dir'{package_dir}'.format(**vars()) )

                shutil.copy(build_path + '/{a}/{a}'.format(**vars()), package_dir)
                if a not in does_not_require_database: dir_util_module.copy_tree(rosetta_dir + '/database', '{package_dir}/database'.format(**vars()), update=False)


            else: raise BenchmarkError('ui_release: ERROR, unsupported os: {platform[os]}!'.format(**vars()))

            release('ui', release_name, package_dir, working_dir, platform, config, release_as_git_repository = False)

            if os.path.isdir(package_dir): shutil.rmtree(package_dir)  # removing package to keep size of database small

        res_code = _S_passed_
        results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : output }
        with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results


def self_release(rosetta_dir, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' self, - do a minimal release procedure to perform a self-test of release system
    '''
    TR = Tracer(True)
    TR('Running self_release test: at working_dir={working_dir!r} with rosetta_dir={rosetta_dir}, platform={platform}, jobs={config[cpu_count]}, memory={config[memory]}GB, hpc_driver={hpc_driver}...'.format( **vars() ) )

    platform_suffix = platform_to_pretty_string(platform)

    release_name = 'self.{platform}.python-{python_version}'.format(platform='.'.join([platform['os']]), python_version=platform['python'].replace('.', '') )

    package_dir = f'{working_dir}/self'
    release_whl  = f'{working_dir}/self.whl'

    os.makedirs(package_dir)
    for i in range(8):
        with open(f'{package_dir}/{i}', 'w') as f: f.write(f'i={i}\n')

    with open(release_whl, 'w') as f: f.write('dummy')

    release('self', release_name, package_dir, working_dir, platform, config, release_as_git_repository = False)

    release('self-wheel', 'self.whl', package_dir=None, working_dir=working_dir, platform=platform, config=config, release_as_git_repository=False, file=release_whl, use_rosetta_versioning=False)

    res_code = _S_passed_
    results = {_StateKey_ : res_code,  _ResultsKey_ : {},  _LogKey_ : 'Done!' }
    with open(working_dir+'/output.json', 'w') as f: json.dump({_ResultsKey_:results[_ResultsKey_], _StateKey_:results[_StateKey_]}, f, sort_keys=True, indent=2)  # makeing sure that results could be serialize in to json, but ommiting logs because they could take too much space

    return results




def py_rosetta4_conda_release(*args, **kwargs): return native_libc_py_rosetta4_conda_release(*args, **kwargs)
#def py_rosetta4_conda_release(*args, **kwargs): return conda_libc_py_rosetta4_conda_release(*args, **kwargs)


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    ''' Run single test.
        Platform is a dict-like object, mandatory fields: {os='Mac', compiler='gcc'}
    '''
    if 'mounts' in config: config['release_root'] = config['mounts']['release_root']

    test = test.replace('PyRosetta4', 'PyRosetta')

    if   test =='source': return rosetta_source_release(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='binary': return rosetta_source_and_binary_release(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug, release_as_git_repository=False)

    elif test =='rosetta.documentation':  return rosetta_documentation(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='PyRosetta.Debug':          return py_rosetta4_release('Debug',          repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.Release':        return py_rosetta4_release('Release',        repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.MinSizeRel':     return py_rosetta4_release('MinSizeRel',     repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.RelWithDebInfo': return py_rosetta4_release('RelWithDebInfo', repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='PyRosetta.conda.Debug':          return py_rosetta4_conda_release('Debug',          repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.conda.Release':        return py_rosetta4_conda_release('Release',        repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.conda.MinSizeRel':     return py_rosetta4_conda_release('MinSizeRel',     repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)
    elif test =='PyRosetta.conda.RelWithDebInfo': return py_rosetta4_conda_release('RelWithDebInfo', repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='PyRosetta.documentation':  return py_rosetta4_documentation('MinSizeRel', repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='ui': return ui_release(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    elif test =='self': return self_release(repository_root, working_dir, platform, config=config, hpc_driver=hpc_driver, verbose=verbose, debug=debug)

    else: raise BenchmarkError('Unknow release test: {}!'.format(test))
