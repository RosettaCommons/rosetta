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

import sys, time, os, re, os.path, subprocess

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

def retrieve_version_information():
    ver = ""
    url = ""
    commit_date = ""

    if subprocess.call("git rev-parse --show-toplevel", shell=True, stdout=subprocess.PIPE) is 0:
        try:
            ver = subprocess.check_output("git rev-parse HEAD", shell=True).strip()
            url = subprocess.check_output("git remote -v | grep fetch | awk '{print $2}' | head -n1", shell=True).strip()
            commit_date = subprocess.check_output("git log %s -1 --format='%%ci'" % ver, shell=True).strip()
        except subprocess.CalledProcessError:
            pass

        if not ver:
            ver = "unknown"
        if not url:
            url = "unknown"
    else:
        # We're probably a release version
        ver = "exported"
        url = "http://www.rosettacommons.org"
        commit_date = ""

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

def generate_version_files():
    '''
    Generates a C++ header file with a summary of the current version(s) of the working copy, if any.
    If this code is not a git repository, the version will be given as "exported".
    Although this is being placed in core/, it doesn't really belong to any subproject.
    There's no good way to know when the version summary will change, either, so we just generate the file every time.
    '''

    version_info = retrieve_version_information()
    # moved to utility/version.hh update_file_if_changed( os.path.normpath("src/devel/svn_version.cc"),  version_cc_template % version_info)
    update_file_if_changed( os.path.normpath("src/python/bindings/src/version.py"), version_py_template % version_info)
    update_file_if_changed( os.path.normpath("src/python/packaged_bindings/src/version.py"), version_py_template % version_info)
    #update_file_if_changed( os.path.normpath("src/python/PyRosetta/src/pyrosetta/version.py"), version_py_template % version_info)

    update_file_if_changed( os.path.normpath("src/utility/version.hh"), version_utility_template % version_info)


def main():
    # Run with timing
    starttime = time.time()
    if '-q' not in sys.argv:
        sys.stdout.write("Running versioning script ... ")
        sys.stdout.flush() # Make sure it gets dumped before running the function.
    generate_version_files()
    if '-q' not in sys.argv:
        sys.stdout.write("Done. (%.1f seconds)\n" % (time.time() - starttime) )

version_py_template = '''\
commit_id = '%(commit_id)s'
commit    = '%(ver)s'
url       = '%(url)s'
date      = '%(commit_date)s'
'''

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
    static std::string commit_id() { return "%(commit_id)s"; }
    static std::string commit()    { return "%(ver)s"; }
    static std::string url()       { return "%(url)s"; }
    static std::string date()      { return "%(commit_date)s"; }
};

} // utility
#endif // INCLUDED_utility_version_hh
'''


version_cc_template = '''// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/svn_version.cc
///
/// @brief
/// @author Ian W. Davis
/// @author Andrew Leaver-Fay
/// @author Sergey Lyskov

/*****************************************************************************
*   This file is automatically generated by the build system.                *
*   DO NOT try modifying it by hand -- your changes will be lost.            *
*   DO NOT commit it to the Git repository.                                  *
*****************************************************************************/

#include <core/svn_version.hh>

namespace devel {

std::string rosetta_svn_version() { return "%(ver)s %(commit_date)s"; }
std::string rosetta_svn_url() { return "%(url)s"; }

class VersionRegistrator
{
public:
	VersionRegistrator() {
		core::set_svn_version_and_url( rosetta_svn_version(), rosetta_svn_url() );
	}
};

// There should only ever be one instance of this class
// so that core::set_svn_version_and_url is called only once
VersionRegistrator vr;

void
register_version_with_core() {
	// oh -- there's nothing in this function.  But
	// forcing devel::init to call this function ensures
	// that the vr variable in this file gets instantiated
}

} // namespace devel
'''

if __name__ == "__main__" or __name__ == "__builtin__":
    os.chdir( os.path.dirname( os.path.abspath(sys.argv[0]) ) ) # where this script is located
    main()
