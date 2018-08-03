#!/usr/bin/env python
# This function expects that the current working directory is the Rosetta root directory.
# If that's ever not true, we need to modify this to take an optional dir name on the cmd line.
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# /// @file   version.py
# ///
# /// @brief
# /// @author Ian W. Davis
# /// @author Andrew Leaver-Fay
# /// @author Sergey Lyskov

from __future__ import print_function

import sys, time, os, re, os.path, subprocess, json, argparse, subprocess, datetime


def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output + errors

        output = output.decode(encoding="utf-8", errors="replace")

        exit_code = p.returncode

        if exit_code  and  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if return_ == 'tuple': return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else: print("\nEncounter error while executing: " + command_line + '\n' + output); sys.exit(1)

    if return_ == 'output': return output
    else: return False


def update_file_if_changed(filename, contents):
    changed = False
    if os.path.exists( filename ):
        f = open(filename);
        try:
            if f.read() != contents:
                changed = True
        finally:
            f.close()
    else:
        changed = True

    if changed:
        f = open(filename, 'w')
        try:
            f.write(contents)
        finally:
            f.close()


def generate_version_information(rosetta_dir, url=None, branch=None, package=None, revision=None, date=None, file_name=None, version=None):
    ''' Generate JSON-like data with Rosetta vesion information. If some of the fields is not supplied either try to autodetect them (if possible) or supply `null` value for it
        Example output:
        {
          "branch"   : "master",                          // Testing server branch
          "revision" : 59804,                             // Testing server revision for this branch
          "version"  : "v2017.42.02-59804+master.a2ef319" // PEP-440 version string
          "package"  : "rosetta.source.master-59804",     // Human redable 'name + version string' for this package
          "url"      : "https://www.rosettacommons.org",  // URL, will be Rosetta or PyRosetta web site address for public releases and `git origin` url for devel releases (ie git@github.com:RosettaCommons/main.git)

          "date": "2017-10-23T20:16:27",
          "week": 35,
          "year": 2017,


          // info for developers

          "source": {
            "binder": "a2df3ed6ea059cf15027ae00a145fc31379272cc",  // Git sha1 for each repository if present
            "demos": "f8e337cddb4089603de5ac2e5d88b00033cdb2f1",
            "documentation": null,
            "main": "f1b74cc7669136f3b0dfcd403fb47726217dac6e",
            "tools": "cc443e94d1d55345e90d96908bfa21575d99143c"
          },

          "git_describe": {                                    // literal output of `git describe --tag` for this git commit and its parsed components
              "commit": "a2ef319",
              "describe": "v2017.42-dev59789-248-ga2ef319",
              "post_revision": 248,
              "revision": 59789,
              "week": 42,
              "year": 2017
          },

        }
    '''
    rosetta_dir = os.path.abspath(rosetta_dir)

    versions = {}
    for repository in 'main tools demos documentation'.split():
        repository_dir = rosetta_dir if repository == 'main' else '{rosetta_dir}/../{repository}'.format(**vars())
        if os.path.isdir(repository_dir):
            res, output = execute('Getting Git commit SHA1 for rosetta {repository}...',  'cd {repository_dir} && git rev-parse HEAD'.format(**vars()), return_='tuple', silent=True)
            versions[repository] = None if res else output[:-1]  # remove \n at the end
        else: versions[repository] = None

    #versions['binder'] = execute('Getting Binder submodule SHA1...', "cd {rosetta_dir} && git ls-tree HEAD source/src/python/PyRosetta/binder | awk '{{print $3}}'".format(**vars()), return_='output')
    for s, p in [('binder', 'source/src/python/PyRosetta/binder'), ('rosetta_scripts_scripts', 'rosetta_scripts_scripts'), ('pyrosetta_scripts', 'pyrosetta_scripts'), ('pybind11', 'source/external/pybind11')]:
        versions[s] = execute('Getting {s} submodule SHA1...'.format(**vars()), "cd {rosetta_dir} && git ls-tree HEAD {p} | awk '{{print $3}}'".format(**vars()), return_='output', silent=True)[:-1]  # remove \n at the end

    if date is None:
        d = execute('Getting date for main {} commit...'.format(versions['main']),  "cd {rosetta_dir} && git log -1 --format='%ci' {versions[main]}".format(**vars()), return_='output', silent=True)
        date = datetime.datetime.strptime(' '.join( d.split()[:2] ), "%Y-%m-%d %H:%M:%S")

    year = date.isocalendar()[0]
    week = date.isocalendar()[1]

    if url is None: url = execute('Getting source origin...',  "cd {rosetta_dir} && git remote -v | grep fetch | awk '{{print $2}}' | head -n1".format(**vars()), return_='output', silent=True)[:-1] # remove \n at the end

    if branch is None: branch = execute('Getting current branch',  'cd {rosetta_dir} && git rev-parse --abbrev-ref HEAD'.format(**vars()), return_='output', silent=True)[:-1] # remove \n at the end

    res, _ = execute('Checking if current Git repositoty is clone of RosettaCommons/main...',  'cd {rosetta_dir} && git cat-file -t e7ed669d70414d073c5477a317a65cea1172daa2'.format(**vars()), return_='tuple', silent=True)
    if res: git_describe = None
    else:
        # Use git-describe --long to always include post-version and sha information
        git_describe_str = execute('Getting `git describe` for current commit...',  "cd {rosetta_dir} && git describe --tags --long --match 'v[0-9]*'".format(**vars()), return_='output', silent=True)[:-1] # remove \n at the end
        describe_match = re.match("v(?P<year>\d+)\.(?P<week>\d+)(-dev(?P<dev_revision>\d+))?-(?P<post_revision>\d+)-g(?P<commit>\w+)", git_describe_str)

        if describe_match:
            git_describe = describe_match.groupdict()
            for int_field in ("year", "week", "dev_revision", "post_revision"):
                if git_describe[int_field] is not None:
                    git_describe[int_field] = int(git_describe[int_field])

            git_describe["describe"] = git_describe_str
        else:
            git_describe = None

    if version is None:
        # 'version' string was not specified, so we assume that we producing local build and Generating PEP-440 version string
        # Expected version string for a release tagged branch is: 2017.42+{branch_names}.{sha1}
        # Expected version string for any other branch is: 2017.42.post.dev+{post_revision.{branch_name}.{sha1}

        mangled_branch = re.sub('[^0-9a-zA-Z.]', '', branch.replace('/', '.'))

        if git_describe and not git_describe["post_revision"]:
            # On a specific tag, so use "release format" tag
            version = '{git_describe[year]}.{git_describe[week]:=02}+{mangled_branch}.{git_describe[commit]}'.format(**vars())
        elif git_describe:
            # On a non-release revision, format a "local" tag
            version = '{git_describe[year]}.{git_describe[week]:=02}.post.dev+{git_describe[post_revision]}.{mangled_branch}.{git_describe[commit]}'.format(**vars())
        else:
            version = 'unknown'

    info = dict(branch = branch,
                revision      = revision,
                package       = package,
                version       = version,
                url           = url,

                date          = date.isoformat(),
                year          = year,
                week          = week,

                source        = versions,
                git_describe  = git_describe,
    )

    if file_name:
        with open(file_name, 'w') as f: json.dump(info, f, sort_keys=True, indent=2)

    return info


def retrieve_version_information():
    json_version_file = None
    release_file = './../.release.json'

    if Options.version:
        print( 'Custom version file was supplied:{}, using it...'.format(Options.version) )
        json_version_file = Options.version
    elif os.path.isfile(release_file):
        print('Release package detected, using rosetta/main/.release.json to acquire version information...')
        json_version_file = release_file

    if json_version_file:
        with open(json_version_file) as f: info = json.load(f)
    else:
        info = generate_version_information(rosetta_dir='./../', file_name='.version.json')

    return info

'''
        ver = ""
        url = ""
        commit_date = ""

        try:
            ver = subprocess.check_output("git rev-parse HEAD", shell=True).strip()
            ver = ver.decode('utf-8', errors="replace") if type(ver) == bytes else ver
            url = subprocess.check_output("git remote -v | grep fetch | awk '{print $2}' | head -n1", shell=True).strip()
            commit_date = subprocess.check_output("git log {} -1 --format='%ci'".format(ver), shell=True).strip()
        except subprocess.CalledProcessError:
            pass

        if not ver:
            ver = "unknown"
        if not url:
            url = "unknown"

        # get_commit_id.sh is not in the standard repository, but is added by PyRosetta?
        # See commit log for  dbbff5655669f41af0dfa7c9421fc89e36b2a227
        commit_id = 'unknown'
        if os.path.isfile('get_commit_id.sh'):
            try:
                res = 0
                output = subprocess.check_output(['./get_commit_id.sh {}'.format(ver)])
            except subprocess.CalledProcessError as err:
                res = err.returncode
                output = err.output

            print('Asked Testing server for commit id, got reply:', repr(output))

            if (res  or  not output  or not output.isdigit() ):
                commit_id = 'failed_to_get_id' # simple validation
            else:
                commit_id = str(int(output))

        if commit_id != 'unknown':
            ver = commit_id + ':' + ver

    return dict( commit_id = commit_id, ver = ver, url = url, commit_date = commit_date)
'''

def generate_version_files():
    '''
    Generates a C++ header file with a summary of the current version(s) of the working copy, if any.
    If this code is not a git repository, the version will be given as "exported".
    Although this is being placed in core/, it doesn't really belong to any subproject.
    There's no good way to know when the version summary will change, either, so we just generate the file every time.
    '''
    info = retrieve_version_information()
    # moved to utility/version.hh update_file_if_changed( os.path.normpath("src/devel/svn_version.cc"),  version_cc_template % version_info)
    # update_file_if_changed( os.path.normpath("src/python/bindings/src/version.py"), version_py_template % version_info)
    # update_file_if_changed( os.path.normpath("src/python/packaged_bindings/src/version.py"), version_py_template % version_info)
    #update_file_if_changed( os.path.normpath("src/python/PyRosetta/src/pyrosetta/version.py"), version_py_template % version_info)

    update_file_if_changed( os.path.normpath("src/utility/version.hh"),
                            version_utility_template % dict(info,
                                                            package = info['package'] if info['package'] else 'devel',
                                                            commit = info['source']['main']) )


# version_py_template = '''\
# package  = '%(package)s'
# version  = '%(version)s'
# revision = '%(revision)s'
# commit   = '%(commit)s'
# url      = '%(url)s'
# date     = '%(date)s'
# '''

version_utility_template = '''\
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/version.hh
///
/// @brief
/// @author Sergey Lyskov

/*****************************************************************************
*   This file is automatically generated by the build system.                *
*   DO NOT try modifying it by hand -- your changes will be lost.            *
*   DO NOT commit it to the Git repository.                                  *
*****************************************************************************/

#ifndef INCLUDED_utility_version_hh
#define INCLUDED_utility_version_hh

#include <string>

namespace utility {

struct Version
{
    static inline std::string package()  { return "%(package)s"; }
    static inline std::string version()  { return "%(version)s"; }
    static inline std::string revision() { return "%(revision)s"; }
    static inline std::string commit()   { return "%(commit)s"; }
    static inline std::string url()      { return "%(url)s"; }
    static inline std::string date()     { return "%(date)s"; }
};

} // utility
#endif // INCLUDED_utility_version_hh
'''

def main(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--quiet', action="store_true", help="Do not output anything during script run.")
    parser.add_argument('--version', help="Supply custom version file. This option override any version detection logic including checks for .release.json file in main/.")

    global Options
    Options = parser.parse_args(args=args[1:])

    starttime = time.time()

    if not Options.quiet: sys.stdout.write("Running versioning script ... "); sys.stdout.flush()

    generate_version_files()

    if not Options.quiet: sys.stdout.write("Done. (%.1f seconds)\n" % (time.time() - starttime) )


if __name__ == "__main__" or __name__ == "__builtin__":
    os.chdir( os.path.dirname( os.path.abspath(sys.argv[0]) ) ) # where this script is located
    main(sys.argv)
