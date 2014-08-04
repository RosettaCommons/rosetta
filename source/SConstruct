# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
#
# Main build file for system.

import os, sys

EnsureSConsVersion(0, 98, 1) # Need version 0.98.1 or above
if( not (2,4) < sys.version_info < (3,0) ):
    print "A version of Python between 2.4 and 2.7 is required to build Rosetta."
    print "Version being used:"
    print sys.version

def print_build(build):
    from tools.build.settings import Settings

    if build.options.log:
        print "Logging information printing for build:", \
              ",".join(build.options.log), "\n"
        if "toplevel" in build.options.log:
            print "Toplevel directory:",
            Settings.write(build.toplevel)
            print "\n"
        if "targets" in build.options.log:
            print "Command-line targets:", ", ".join(COMMAND_LINE_TARGETS) or "none"
            print "Build targets:", ", ".join(BUILD_TARGETS) or "defaults"
            print ""
        if "platform" in build.options.log:
            print "Platform:",
            Settings.write(build.platform)
            print "\n"
            if 'unit_test_platform_only' in build.targets:
                sys.exit()
        if "projects" in build.options.log:
            print "Projects:",
            Settings.write(build.projects, 1)
            print "\n"
        if "options" in build.options.log:
            print "Command-line options:", \
                  ", ".join( [ "%s = %s" % item for item in ARGUMENTS.items() ] )
            print ""
            print "Requested options:",
            Settings.write(build.options_requested, 1)
            print "\n\nActual options:",
            Settings.write(build.options, 1)
            print "\n"
        if "settings" in build.options.log:
            print "Settings ids:",
            Settings.write([ s.id for s in build.settings ], 1)
            print "\n\nSettings:",
            Settings.write(build.settings, 1)
            print "\n"
        if "environment" in build.options.log:
            print "Environment:"
            env_vars = []
            env_dict = build.environment.Dictionary()
            env_keys = env_dict.keys()
            env_keys.sort()
            for key in env_keys:
                if env_dict.has_key(key):
                    env_vars += [ key + " = " + str(env_dict[key]) ]
            Settings.write(env_vars, 1)
        print "-" * 72

help = """
This is the Rosetta build system (using SCons).
For information on using it properly see the file tools/build/README.
If you are a developer also take a look at the Rosetta wiki under 'Build System.'

The basic form of a run is:

    scons [flags] [<options>] [<targets>]

[flags] are command-line options for the scons program (see more
by using -H).

<options> are key=value pairs from the following sets in any order.

    cxx:        The C++ compiler
    cxx_ver:    The version of the C++ compiler
    os:         The operating system
    os_ver:     The version of the operating system
    arch:       The processor architecture
    arch_size:  The bit-size of the processor
    mode:       debug, release, release_debug, profile or coverage build
    cat:        Whether the build is for sources, tests, and/or docs.
    extras:     Additional features of the build

<targets> are lists of projects or directories you want to build.

Examples:

    scons
        Build the default projects with default settings (debug mode, shared libs)

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

    scons bin pilot_apps_all
        Build all pilot_apps listed in src/pilot_apps.src.settings.all
        The pilot_apps and devel projects must be listed in the either
        tools/build/user.settings or in tools/build/site.settings by
        adding  the following line into the "appends" sections

        "projects" : { "src" : ["pilot_apps","devel"], },

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

    scons extras=apbs

    scons extras=hdf5
        Enable support for HDF5 datastores.
"""

def main():
    DEBUG = True
    try:
        build = SConscript("tools/build/setup.py")
        build.toplevel = os.getcwd()

        print_build(build)

        # Put all .sconsign files (used to determine if a targets
        # has changed) into the build directory
        ### SConsignFile("build/signatures")
        SConsignFile()

        Help(help)

        # Run the SConscript to choose which subdirectories to invoke
        SConscript("SConscript", "build")

    except StandardError, ex:
        if not DEBUG:
            print "Error:", ex
        else:
            # Print a stack trace
            sys.excepthook(*sys.exc_info())

main()
