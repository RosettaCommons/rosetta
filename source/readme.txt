(c) Copyright Rosetta Commons Member Institutions.
(c) This file is part of the Rosetta software suite and is made available under license.
(c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
(c) For more information, see http://www.rosettacommons.org. Questions about this can be
(c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

<< Building >>

In order to build a debug version of the mini-rosetta executables, simply run
scons like this:

scons bin

In order to build release executables, add the flag 'mode=release' like this:

scons bin mode=release

To display more rosetta specific build options call
scons -h
To display general scons build options call
scons -H

Instructions for native builds on OSX using XCode and Windows using Visual
Studio to follow soon.

<< Documentation >>

Mini Rosetta wiki home page: https://wiki.rosettacommons.org/index.php/Projects:miniRosetta
Nightly build of Mini documentation: https://www.rosettacommons.org/internal/doc/mini/index.html
libRosetta Docs: https://wiki.rosettacommons.org/index.php/Projects:Program_Options

<< HowTos >>

Option system: https://wiki.rosettacommons.org/index.php/Projects:mini:Options_HowTo
Unit tests: https://wiki.rosettacommons.org/index.php/Tools:Unit_Tests

<< Additional Build Environment Setup >>

Automatic location of other compilers (assuming they are already in
your path) such as Intel C/C++ may be enabled by uncommenting the
"program_path" line in 'tools/build/user.settings'.

A user can restrict compilation of the the devel and pilot_apps. On issuing
the call

scons bin my
or
scons bin my_pilot_apps

SCons will read from src/devel.src.settings and src/pilot_apps.src.settings
rather than src/devel.src.settings.all and src/pilot_apps.src.settings.all
This cuase SCons to build only the sources listed in src/devel.src.settings
and src/pilot_apps.src.settings with the needed dependencies.

<< Common build calls that may be useful >>

    scons
        Build the default projects with default settings (debug mode, shared
libs)

    scons <project>
        Build the target <project> with default settings

    scons <project>/<subdirectory>
        Build only the sources of <project> in <subdirectory>

    scons <project>/<path/<objectfile>
        Build only <objectfile>

    scons bin
        Build and install executables in bin/ directory

    scons -D #bin
        Build and install executables in bin/ directory if current
        working directory is a sub-directory mini. -D options tells
        scons to iteratively search towards the root for SConstruct
        file. The # sign is an alias for the top build directory.  

    scons bin
    scons bin pilot_apps_all
        Build all pilot_apps listed in src/pilot_apps.src.settings.all
        and sources in src/devel.src.settings as well as the core, numeric
        and utility libraries

    scons bin my
    scons bin my_pilot_apps
        Build restricted set of pilot_apps and devel sources listed in 
        src/pilot_apps.src.settings.my and src/devel.src.settings.my

    scons bin/exec
        Build and install a particular executable in the bin directory
        e.g
        scons bin/benchmark.linuxgccdebug
        scons mode=release bin/benchmark.linuxgccrelease

    scons mode=release
        Build in release mode (~10x faster executable)

    scons extras=static
        Static linking instead of shared libraries (more portable)

    scons cat=test
    python test/run.py
        Build and run unit tests. (Note the sources must be built first.)

    scons -j3
        Parallelize build into 3 threads (faster on multiproc. machine)

    python external/scons-local/scons.py
        Use the version of scons that is distributed with mini
        (Hint: use if scons is not installed on system)
        
+5

