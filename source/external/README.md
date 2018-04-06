EXTERNAL LIBRARIES FOR ROSETTA
==============================

These are libraries which are used by Rosetta, but are not written by RosettaCommons members.

See the individual library directories for details on authorship and licensing of each library.


Adding new external libraries
-----------------------------

0) Licensing

Before adding any external library to Rosetta, please be careful that the license of the
external library is compatible with the Rosetta licence. DO NOT add any library with an
incompatible license to the main Rosetta distribution.

1) Header-only libraries

The ideal external library is one which is header only, and does not require any compilation.
To add a new header library, all you need to do is add it to the main/source/external/ directory.
By default, both the main/source/external/ and main/source/external/include/ directories
are in the include search path.


2) Libraries which require compilation.

To add a library which requires compilation, there are two options for specifying which files
are to be compiled.

In both cases, you will need to add your library name to the "external" section of the
main/source/projects.settings file to tell the build system it exists.


2a) Simple compilation

If you can use Rosetta-like build settings for your library, you can have a
main/source/external/<projectname>.external.settings file for your library.

This will use the same settings as Rosetta, but without any of the warning settings.

Additionally, the *.external.settings file can contain a "defines" list, for enabling
particular preprocessor defines while compiling the external library specifically.

There is also a "only_with_extras" list, which will compile the external library only
when one of those particular extras builds are enabled. Please note that when library specify
`only_with_extras` it will still be build and link for default-build as well, but for
default-build it will be build as library with an empty number of source files.

2b) Custom compilation

If you need special compilation settings, you can add a SConscript.external.<projectname> file to
the main/source/external/ directory. If present, this script control the Scons settings.

You will also need to change the build process in main/source/cmake/build/{external.cmake,external.static.cmake}
to also implement your custom build settings for the CMAKE build system.
