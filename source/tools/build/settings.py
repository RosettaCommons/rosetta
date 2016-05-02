# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
#
# Settings are the fundamental building blocks of the build system.
# They are essentially just dictionaries, but with restricted keys
# and using object syntax for ease of reading (and fewer mysterious
# bugs from typos.)  Their most important feature is how they are
# converted into SCons.Environmen variables.
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import copy, re, string, sys
import SCons # needed for version checking

class Settings(dict):
    """A dictionary which allows its internal items to be accessed
as fields of an object.  The object is limited to just those fields
to catch mistakes (via the Python >= 2.2 "slots" mechanism.
"""


    # This is a small hack to allow us to assign to self.rest
    # in this class.  __slots__ will be different for each
    # subclass.
    __slots__ = [ "rest" ]

    def __init__(self, init = None, **settings):
        dict.__init__(self)
        # Default initialize all fields
        for field in self.__slots__:
            self[field] = copy.copy(init)
        self.rest = {}
        # Add all settings which are known to the object's fields.
        # Add all remaining settings to an external dictionary.
        for name, value in settings.items():
            if name in self.__slots__:
                self[name] = value
            else:
                self.rest[name] = value

    def __getattr__(self, key):
        if key in self.__slots__:
            try:
                return self[key]
            except KeyError:
                raise AttributeError, "No such attribute " + key
        elif key in ("rest",):
            return dict.__getattr__(self, key)
        else:
            raise AttributeError, "No such attribute " + key

    def __setattr__(self, key, value):
        if key in self.__slots__:
            try:
                self[key] = value
            except KeyError:
                raise AttributeError, "No such attribute " + key
        elif key in ("rest",):
            dict.__setattr__(self, key, value)
        else:
            raise AttributeError, "No such attribute " + key

    def keys(self):
        """The order of keys in a settings object is the same
as the order of its slots.
"""
        return self.__slots__

    def iterkeys(self):
        """The order of keys in a settings object is the same
as the order of its slots.
"""
        return iter(self.__slots__)

    def load(class_, file, root = None):
        """Read in the settings from a file.  The file format
is defined to be abstract and accessible ONLY through this
method, which will return the values as a map.  In practice
the settings files are very deliberately Python data (maps,
lists, strings, and numbers) and nothing else.  It is almost
but not quite JSON (www.json.org).
"""
        settings = {}
        execfile(file, settings)
        if settings.has_key("__builtins__"):
            del settings["__builtins__"]
        if root:
            settings = settings[root]
        return settings

    load = classmethod(load)

    def write(class_, settings, indent = 0, file = sys.stdout):
        """Pretty-print the items in a settings object (or in
fact any Python basic type.
(This really should be called "print" not "write" but there's
obviously a keyword conflict.)
"""
        # Settings are also dictionaries
        if isinstance(settings, dict):
            if len(settings) == 0:
                print >>file, "{}",
            else:
                keys = settings.keys()
                for key in keys:
                    value = settings[key]
#                   # Don't print empty items
#                   if not (\
#                       type(value) is type(None) or \
#                       type(value) is str and len(value) == 0
#                   ):
                    if True:
                        print >>file, "\n" + (indent * 2) * " " + key + ":",
                        class_.write(value, indent + 1, file)
        elif type(settings) in (list, tuple):
            if len(settings) == 0:
                print >>file, "[]",
            for item in settings:
                print >>file, "\n" + (indent * 2) * " " + "-", # str(count) + ":",
                class_.write(item, indent + 1, file)
        elif type(settings) is type(None):
            print >>file, "<NONE>",
        elif type(settings) is str and len(settings) == 0:
            print >>file, "\"\"",
        else:
            print >>file, str(settings),

    write = classmethod(write)


class BuildOptionsSupported (Settings):
    """Table of supported options.  This treats versions and bit-sizes as
subsidiary to the settings listed here.
"""

    __slots__ = [
        "cxx", "os", "arch", "mode", "cat", "extras", "log",
    ]

    def __init__(self, options):
        Settings.__init__(self, **options)


class BuildOptions (Settings):
    """The information needed to select the appropriate build targets.
Modifiable either on the command-line or in a .settings file.
Build options are essentially descriptions of the platform: the compiler,
OS, architecture, build mode, and build cat that will control what options
to pass to the compiler.
    """

    __slots__ = [
        "cxx", "cxx_ver", "os", "os_ver", "arch", "arch_size",
        "mode", "cat", "extras", "log",
    ]

    def __init__(self, options):
        Settings.__init__(self, **options)


class BuildFlags (Settings):
    """Types of compiler flags, which it is often helpful to break
up into categories for explanatory purposes.
"""

    __slots__ = [
        "cc", "cxx", "compile", "link", "shlink", "mode", "warn",
    ]

    def __init__(self, **flags):
        Settings.__init__(self, init = [], **flags)
        for name, value in self.items():
            assert type(value) is list, "Parameters of BuildFlags must be lists"


class BuildSettings (Settings):
    """Settings are flags and options passed to the compiler or the
SCons environment to control the build.  They are specified in .settings
files and are combined algorithmically to form the final environment
used by the build.
"""

    __slots__ = [
        "projects", "program_path",
        "include_path", "includes",
        "library_path", "libraries",
        "cc", "cxx", "flags", "defines", "version", "shlinkcom"
    ]

    def __init__(self, settings):
        Settings.__init__(self, "", **settings)
        if type(self.flags) is dict:
            self.flags = BuildFlags(**self.flags)

    def parse_major_minor_revision(self, version_string):
        """Split a version string into major, minor and (optionally)
        revision parts.

        This is complicated by the fact that a version string can be
        something like 3.2b1."""
        version = string.split(string.split(version_string, ' ')[0], '.')
        v_major = int(version[0])
        v_minor = int(re.match('\d+', version[1]).group())
        if len(version) >= 3:
            v_revision = int(re.match('\d+', version[2]).group())
        else:
            v_revision = 0
        return v_major, v_minor, v_revision

    def symbols(self):
        """Generate a dictionary with the symbol names SCons understands.
        """

        # we need to grab the scons version as it will dictate how
        # flags are put together
        scons_version = self.parse_major_minor_revision( SCons.__version__ ) # tuple of ( major, minor, revision )

        symbols = {}
        if self.flags:
            flags = self.flags.compile + self.flags.warn + self.flags.mode
            if self.includes:
                flags += [ "include %s" % (flag) for flag in self.includes ]
            flags = [ "-%s" % (flag) for flag in flags ]

            # store flags for C/C++ compiler
            compiler_flags = flags # flags passed to both C and C++ compilers
            c_flags = [ "-%s" % (flag) for flag in self.flags.cc ] # flags passed to just C compiler
            cxx_flags = [ "-%s" % (flag) for flag in self.flags.cxx ] # flags passed to just C++ compiler

            link_flags = [ "-%s" % (flag) for flag in self.flags.link ]
            link_flags = [ flag.replace("-$", "$") for flag in link_flags ]
            shlink_flags = [ "-%s" % (flag) for flag in self.flags.shlink ]
            shlink_flags = [ flag.replace("-$", "$") for flag in shlink_flags ]
        elif self.includes:
            flags = [ "-include %s" % (flag) for flag in self.includes ]
            compiler_flags = flags
            c_flags = []
            cxx_flags = []
            link_flags = []
            shlink_flags = []
        else:
            compiler_flags = []
            c_flags = []
            cxx_flags = []
            link_flags = []
            shlink_flags = []
        if self.cc:
            symbols["CC"] = self.cc # + ".".join(self.version)
        if self.cxx:
            symbols["CXX"] = self.cxx # + ".".join(self.version)

        # symbol setup here must be broken down by scons version
        # due to scons bugs #846 and #1971
        if scons_version < ( 0, 96, 94 ): # bug and 846 and #1971
            if compiler_flags: # setup flags for both C and C++ compiler
                symbols["CCFLAGS"] = " " + " ".join(compiler_flags) + " " + " ".join(c_flags) # CFLAGS doesn't exist
                symbols["CXXFLAGS"] = " " + " ".join(compiler_flags) # CCFLAGS not passed properly with CXXFLAGS
            if cxx_flags: # setup flags for just C++ compiler
                if "CXXFLAGS" not in symbols:
                    symbols["CXXFLAGS"] = ""
                symbols["CXXFLAGS"] = symbols["CXXFLAGS"] + " " + " ".join(cxx_flags)
        elif scons_version < ( 0, 98, 1 ): # bug #1971
            if compiler_flags: # setup flags for both C and C++ compiler
                symbols["CCFLAGS"] = " " + " ".join(compiler_flags)
                symbols["CXXFLAGS"] = " " + " ".join(compiler_flags) # CCFLAGS not passed properly with CXXFLAGS
            if c_flags: # setup flags for just C compiler
                symbols["CFLAGS"] = " " + " ".join(c_flags)
            if cxx_flags: # setup flags for just C++ compiler
                if "CXXFLAGS" not in symbols:
                    symbols["CXXFLAGS"] = ""
                symbols["CXXFLAGS"] = symbols["CXXFLAGS"] + " " + " ".join(cxx_flags)
        else:
            if compiler_flags: # setup flags for both C and C++ compiler
                symbols["CCFLAGS"] = " " + " ".join(compiler_flags)
            if c_flags: # setup flags for just C compiler
                symbols["CFLAGS"] = " " + " ".join(c_flags)
            if cxx_flags: # setup flags for just C++ compiler
                symbols["CXXFLAGS"] = " " + " ".join(cxx_flags)

        if self.defines:
            symbols["CPPFLAGS"] = [ "-D%s" % d for d in self.defines ] # " " + " ".join(self.defines)
        if link_flags:
            symbols["LINKFLAGS"] = " " + " ".join(link_flags)
        if shlink_flags:
            symbols["SHLINKFLAGS"] = " " + " ".join(shlink_flags)
        if self.shlinkcom:
            symbols["SHLINKCOM"] = self.shlinkcom
        if self.include_path:
            symbols["CPPPATH"] = self.include_path
        # self.includes is handled with self.flags
        if self.library_path:
            symbols["LIBPATH"] = self.library_path
        if self.libraries:
            symbols["LIBS"] = self.libraries
        if self.program_path:
            symbols["ENV"] = {}
            symbols["ENV"]["PATH"] = self.program_path
        symbols.update(self.rest)
        return symbols



class BuildSettingsCombined (Settings):
    """The settings which are merged to create the final environment SCons
uses for it's build.  Combined settings are configured from a settings file.
They store four sets of BuildSettings, of which the first is by far the
most common:
  - appends, which append their settings to the current environment.
  - prepends, which prepend their settings to the current environment.
  - overrides, which replace any existing settings of the same name in the
    current environment.
  - removes, which delete any existing settings which match their value.
"""

    __slots__ = [ "id", "appends", "prepends", "overrides", "removes" ]

    def __init__(self, id, settings):
        Settings.__init__(self, **settings)
        self.id = id
        for name in self.__slots__[1:]:
            value = getattr(self, name)
            assert type(value) is dict or type(None)
            if type(value) is dict:
                setattr(self, name, BuildSettings(value))


class ProjectSettings(Settings):
    """Settings for an individual project.  Added to build settings
by SConscripts.
    """

    __slots__ = [
        "name", "title", "version",
        "subprojects",
        "sources", "extras",
        "include_path", "includes",
        "library_path", "libraries",
        "defines", "testinputfiles",
        "ccflags","cflags","cxxflags","link_flags",
    ]

    def __init__(self, name, settings = None):
        if not settings: settings = {}
        Settings.__init__(self, init = [], **settings)
        self.name = name

    def _get_packages(self):
        packages = self.sources.keys()
        packages.sort()
        # packages = [ self.name + "/" + package for package in packages if package ]
        return packages

    packages = property(_get_packages)

    def _get_modules(self):
        modules = []
        for package in self.packages:
            if package: slash = "/"
            else:       slash = ""
            modules += [ "%s%s%s" % (package, slash, module)
                         for module in self.sources[package] ]

        return modules

    modules = property(_get_modules)

    def modules_for(self, package):
        return self.sources[package]

    def symbols(self):
        symbols = {}
        if self.defines:
            symbols["CPPFLAGS"] = self.defines # " ".join(self.defines) + " "
        if self.include_path:
            symbols["CPPPATH"] = self.include_path
        if self.includes:
            symbols["CCFLAGS"] = self.includes
            symbols["CXXFLAGS"] = self.includes
        if self.library_path:
            symbols["LIBPATH"] = self.library_path
        if self.libraries:
            symbols["LIBS"] = self.libraries
        symbols.update(self.rest)
        return symbols


class BuildState (Settings):
    """Global state for a build.
"""

    __slots__ = [
        "environment",
        "toplevel", "projects",
        "platform", "platform_includes",
        "options", "options_requested",
        "settings", "targets", "all_libraries"
    ]

    def __init__(self, state = None):
        if not state: state = {}
        Settings.__init__(self, **state)


