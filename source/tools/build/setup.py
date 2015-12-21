# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
#
# Main build file for system.
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os, os.path, sys

from tools.build import utility, setup_platforms
from tools.build.settings import Settings, BuildState, \
    BuildOptions, BuildOptionsSupported, BuildSettings, BuildSettingsCombined


def copy_template_files():
    '''Copy *.template to * (if not exists) to help new users get set up'''
    for filename in ["user.options", "user.settings", "#src/pilot_apps.src.settings", "#src/devel.src.settings"]:
        src = File(filename+".template")
        dst = File(filename)
        if not os.path.exists( str(dst) ):
            Execute( Copy(str(dst), str(src)) ) # yes, SCons lists them in reverse order

def setup_toplevel():
    toplevel = os.path.normpath(os.getcwd() + "/../..")
    return toplevel

def setup_build_options():
    """Select the build options specified, either on the command-line
or in .options files, which determine the type of the build, what
directory it is built to, and what settings it ultimately uses.
    """

    # Import supported options
    supported = Settings.load("options.settings", "options")
    supported = BuildOptionsSupported(supported)

    # Import default options
    defaults = Settings.load("basic.options")
    defaults = BuildOptions(defaults)

    # Incorporate options in site and user options
    option_files = [ "%s.options" % (name) for name in ("site", "user") ]
    options = Variables(option_files)

    # Read in the options from the .settings files (1st pass)
    # This has to be done in two passes, because some legal values
    # (version numbers) are dependent on the other values
    options.AddVariables(
        # Options for the build.  Modifiable on the command-line.
        EnumVariable("cxx", "Select the C++ compiler to build with",
                     defaults.cxx, supported.cxx),
        EnumVariable("os", "Select the operating system to build for",
                     defaults.os, supported.os),
        EnumVariable("arch", "Select the processor architecture to build for",
                     defaults.arch, supported.arch),
        EnumVariable("mode", "Select the optimization mode",
                     defaults.mode, supported.mode),
        ListVariable("cat", "Select the build category",
                     defaults.cat, supported.cat),
        ListVariable("extras", "Select any extras",
                     defaults.extras, supported.extras),
        ListVariable("log", "Log debugging output",
                     defaults.log, supported.log),
    )

    # Incorporate dependent values from .settings files (2nd pass)
    # Creating an Environment is necessary because it's not possible
    # to extract values from an Options object directly.  This
    # environment will be discarded.
    env  = Environment(options = options)
    cxx  = env["cxx"]
    os   = env["os"]
    arch = env["arch"]
    options.AddVariables(
        EnumVariable("cxx_ver", "Select the C++ compiler version",
                     defaults.cxx_ver, supported.cxx[cxx]),
        EnumVariable("os_ver", "Select the operating system version",
                     defaults.os_ver, supported.os[os]),
        EnumVariable("arch_size", "Select the architecture bit-size",
                     defaults.arch_size, supported.arch[arch]),
        )

    # Create the environment with the complete options.
    # This is mostly for validation of the options: the values
    # of the options will be extracted into 'requested' and 'actual'
    # BuildOptions, and this environment will be discarded.
    env = Environment(options = options)

    # Generate a build options object for the options requested.
    # This is used primarily for error reporting.
    # (There must be some dict function for extracting a subset of a map...)
    requested = BuildOptions(
        utility.map_subset(env.Dictionary(), defaults.keys())
    )

    # Generate a build options object for the actual options used.
    # This is used to select the build targets.
    actual = BuildOptions(requested)

    # Verify/validate that all of the options are values that are
    # meaningful target id as specified in targets.settings.

    # The C++ compiler
    actual.cxx = setup_platforms.select_compiler(supported, requested.cxx)
    actual.cxx_ver = setup_platforms.select_compiler_version(
        supported, actual.cxx, requested.cxx_ver
    )

    # The operating system
    actual.os = setup_platforms.select_os(supported, requested.os)
    actual.os_ver = setup_platforms.select_os_version(
        supported, actual.os, requested.os_ver
    )
    # The platform architecture
    actual.arch = setup_platforms.select_arch(
        supported, actual.os, requested.arch
    )
    actual.arch_size = setup_platforms.select_arch_size(
        supported, actual.os, actual.arch, requested.arch_size
    )
    # The build cat
    actual.cat = setup_platforms.select_cat(supported, requested.cat)

    # Nothing special is currently done to select mode, kind or extras
    if len(actual.extras) == 0:
        actual.extras.append("default")
    else :
        actual.extras.data.sort()

    return requested, actual



def setup_build_settings(options):
    """Select the build settings needed to generate the combined settings.
The order of merging of settings is set programmatically rather than explicitly,
because explicit enumeration of all possible build types becomes rapidly
combinatorically unmanageable.
"""

    # Import default settings
    supported = Settings.load("basic.settings", "settings")

    # Create ids from a combination of options
    possible = [ "base" ]

    # Build options which are strings can be turned into ids
    # fairly systematically.  The exceptions come afterward.
    # Build options which are lists (i.e. cat and extras) need to be
    # incorporated separately for every item in the option.

    # OS and compiler base lines
    for keys in [
        ("os",), ("os", "os_ver"),
        ("cxx",), ("cxx", "cxx_ver"),
    ]:
        values = [ getattr(options, name) for name in keys ]
        possible += [ ", ".join(values) ]

    # Compilers and cat
    for kind in options.cat:
        # ("cxx", "kind")
        possible += [ ", ".join((options.cxx, kind)) ]

    # Compiler, OS, architecture and mode combinations
    for keys in [
        ("cxx", "os"), ("cxx", "cxx_ver", "os"),
        ("cxx", "os", "os_ver"), ("cxx", "cxx_ver", "os", "os_ver"),
        ("cxx", "arch"), ("cxx", "cxx_ver", "arch"),
        ("cxx", "arch", "arch_size"), ("cxx", "cxx_ver", "arch", "arch_size"),
        ("cxx", "os", "arch", "arch_size"),
        ("cxx", "mode"), ("cxx", "cxx_ver", "mode"),
        ("cxx", "os", "mode"),  ("cxx", "cxx_ver", "os", "mode"),
        ("cxx", "arch", "mode"),  ("cxx", "cxx_ver", "arch", "mode"),
        ("cxx", "os", "arch"), ("cxx", "os", "arch", "mode" ),
    ]:
        values = [ getattr(options, name) for name in keys ]
        possible += [ ", ".join(values) ]

    # Compilers and extras
    for extra in options.extras:
        # ("cxx", "extra")
        possible += [ ", ".join((options.cxx, extra)) ]
        # ("cxx", "os", "extra")
        possible += [ ", ".join((options.cxx, options.os, extra)) ]
        # ("cxx", "arch", "extra")
        possible += [ ", ".join((options.cxx, options.arch, extra)) ]
        # ("cxx", "cxx_ver", "extra")
        possible += [ ", ".join((options.cxx, options.cxx_ver, extra)) ]


    # Setting up paths for Windows and MinGW by importing PATH from the environment
    # Windows builds have no canonical system location for libraries.
    if options.os in "windows" and "mingw" in options.extras:
        sys_path = os.environ["PATH"]
        windows_or_mingw = {
            "windows, mingw" : {
                "overrides" : {
                    "program_path" : sys_path,
                },
            },
        }
        supported.update(windows_or_mingw)
        possible += [ "windows, mingw" ]

    # Site settings
    if (os.path.exists("site.settings")):
        site = Settings.load("site.settings", "settings")
        supported.update(site)
    possible += [ "site" ]

    # User settings
    if (os.path.exists("user.settings")):
        user = Settings.load("user.settings", "settings")
        supported.update(user)
    possible += [ "user" ]

    # Actual ids are those possible ids which are actually supported
    # by the settings in the settings file.
    actual = []
    for id in possible:
        # It's not an error not to have defined a given id
        # We only incorporate the settings of those ids we happen to find.
        # This could lead to subtle errors, because there's no easy way
        # to recognize the absence of a given id.
        if supported.has_key(id) and supported[id]:
            actual += [ id ]

    # Create the build settings which match all actual build ids
    settings = []
    for id in actual:
        settings += [ BuildSettingsCombined(id, supported[id]) ]
        # Debug messages for checking evaluation order
        #print str(id) + ":"
        #for key, value in supported[id].iteritems():
        #    print "  " + key + ":"
        #    for key, value in value.iteritems():
        #        print "    " + key + ": " + str(value)

    return settings


def setup_platform_path(options):
    """Generate the platform component of the build path.
These are of the form:
    <mode>/<os>/<os_ver>/<arch_size>/<arch>/<cxx>/<cxx_ver>[/<extras>]
e.g:
    debug/linux[-2.6]/x86[-32]/gcc-4.0.0[/boinc]
    """
    if options.cxx_ver == "*":
        cxx_ver = ""
    else:
        cxx_ver = options.cxx_ver
    path = [ options.mode,
             options.os, options.os_ver,
             options.arch_size, options.arch,
             options.cxx, cxx_ver,
             "-".join(options.extras)
           ]
    path = os.path.join(*path).rstrip(setup_platforms.get_os_pathsep())
    return path


def setup_platform_includes(options):
    """Generate the path to the platform dependent includes.
These are of the form:
    <os>/<os_ver>/<arch_size>/<arch>/<compiler>
"""
    includes = [ options.os, \
#                 options.os_ver, \        #SGM Simplified the platform tree
                 options.arch_size,
#                 options.arch, \          #SGM Simplified the platform tree
                 options.cxx ]
    if options.cxx_ver != "*":
        includes += [ options.cxx_ver ]
    # includes = os.path.join(*includes)
    return includes


def get_projects():
    """Get a structure of all the allowed projects."""

    # Import project names from the settings files
    allowed = Settings.load("../../projects.settings", "projects")

    #print DEBUG
    #print allowed

    #The following block allows users to add, remove, and override which projects are built by default
    #on the site and user level by including the "project" keyword in the user.settings and site.settings file
    modifications = [ "site", "user"]
    settings_exist = 0
    for modification in modifications :
        import os
        if os.path.exists( modification + ".settings") :
            site = Settings.load( modification + ".settings", "settings");

            if site[modification]["prepends"].has_key("projects") :
                for cat in site[modification]["prepends"]["projects"] :
                    for item in site[modification]["prepends"]["projects"][cat] :
                        allowed[cat].insert(0, item )
            if site[modification]["appends"].has_key("projects") :
                for cat in site[modification]["appends"]["projects"] :
                    for item in site[modification]["appends"]["projects"][cat] :
                        allowed[cat].append(item)
            if site[modification]["overrides"].has_key("projects") :
                for cat in site[modification]["overrides"]["projects"] :
                    allowed[cat] = site[modification]["overrides"]["projects"][cat]
            if site[modification]["removes"].has_key("projects") :
                for cat in site[modification]["removes"]["projects"] :
                    for project in allowed[cat] :
                        for removed_project in site[modification]["removes"]["projects"][cat] :
                            if removed_project == project :
                                allowed[cat].remove(project)

    #Print allowed projects
    #print allowed

    return allowed

def setup_projects(targets):
    """The projects are not specified as keyword arguments but
as unadorned ones, which SCons stores as COMMAND_LINE_TARGETS.
Check to see that all command line targets are ones which are
legitimately in the project settings, and apply those that are
as the final project list.
    """

    actual = []
    invalid = []

    allowed = get_projects()

    # Keep those project names which are found in the path
    # elements of the targets.
    # XXX: This is not really reliable.  If a project name is
    # part of the path of another project it will be falsely
    # accepted.
    actual = []
    for project in allowed :
        for target in targets:
            elements = target.split("/")
            if project in elements:
                actual += [ project ]

    if not actual:
        actual = allowed

    return actual


def setup_environment(settings):
    """Setup the SCons environment.
Do this by combining all the symbols in each build setting appropriately
into the fresh environment.  When this is complete the environment is
what SCons will use to build the system.
"""
    # CXXFLAGS have to be cleared because SCons dumps all the CCFLAGS
    # into them.
    env = Environment()
    for setting in settings:
        if setting.appends:
            symbols = setting.appends.symbols()
            if symbols.has_key("ENV"):
                for key, value in symbols["ENV"].items():
                    env["ENV"][key] += value
                del symbols["ENV"]
            env.Append(**symbols)
        if setting.prepends:
            symbols = setting.prepends.symbols()
            if symbols.has_key("ENV"):
                for key, value in symbols["ENV"].items():
                    if key in env["ENV"]:
                        env["ENV"][key] = value + env["ENV"][key]
                    else:
                        env["ENV"][key] = value
                del symbols["ENV"]
            env.Prepend(**symbols)
        if setting.overrides:
            env.Replace(**setting.overrides.symbols())
        if setting.removes:
            removes = setting.removes.symbols()
            for name, values in removes.items():
                existing = env.Dictionary()[name]
                if isinstance( existing, list ): # list of values
                     for value in values:
                          if value in existing:
                               existing.remove(value)
                else: # everything else
                     # A little extra work here because flag values in
                     # env.Dictionary() are being stored as strings not lists.
                     # See BuildSettings.symbols() for further comments.
                     existing = str(env.Dictionary()[name])
                     existing = existing.split()
                     values = values.split()
                     for value in values:
                         if value in existing:
                             existing.remove(value)
                     env.Dictionary()[name] = " ".join(existing)

    # Use pipes instead of shell commands so that SCons doesn't swallow the output
    # from running a compiler or tool.
    env["PIPE_BUILD"] = True
    env["PSTDOUT"] = sys.stdout
    env["PSTDERR"] = sys.stderr
    return env


def setup():
    """Combine all setup functions to initialize build state.
"""
    copy_template_files()
    build = BuildState()
    # Store the toplevel directory
    # (SCons grants access using the '#' token but doesn't store the path)
    # build.toplevel = setup_toplevel()
    build.options_requested, build.options = setup_build_options()
    build.settings = setup_build_settings(build.options)
    build.environment = setup_environment(build.settings)
    build.platform = setup_platform_path(build.options)
    build.platform_includes = setup_platform_includes(build.options)
    build.all_libraries = get_projects()
    build.projects = setup_projects(COMMAND_LINE_TARGETS)
    build.targets = BUILD_TARGETS
    return build

build = setup()
Return("build")
