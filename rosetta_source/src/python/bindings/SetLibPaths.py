#!/usr/bin/env python
"""This set of functions allows the recursive setting of the dynamic library
names so that they can find one another regardless of the absolute location
in the file system. (Mac OS X only)

Colin A. Smith <colin.smith@ucsf.edu>"""

import subprocess
import os.path
import re
import sys

VERBOSE = False
ROSETTA_LIBS = ["libObjexxFCL.dylib", "libutility.dylib", "libnumeric.dylib",
                "libcore.dylib", "libprotocols.dylib", "libdevel.dylib"]

def getLibraryNames(filepath):
	"""get library names using otool -L"""

	output = subprocess.Popen(["/usr/bin/otool", "-L", filepath], stdout=subprocess.PIPE).communicate()[0]

	return [line.strip().split(" (")[0] for line in output.split("\n")[1:-1]]

def setLibraryNames(filepath, dict = {}, id = None):
	"""set library names using install_name_tool"""

	install_name_tool_args = ["/usr/bin/install_name_tool"]

	for key, value in dict.iteritems():
		install_name_tool_args += ["-change", key, value]

	if id:
		install_name_tool_args += ["-id", id]

	install_name_tool_args += [filepath]

	if VERBOSE:
		print " ".join(install_name_tool_args)
	subprocess.Popen(install_name_tool_args).communicate()

def setLibraryPath(filepath, libfilenames, libpath = "", loader_path = False):
	"""set a subset of the dependent library names in a dynamic library to
	the given path"""

	libnames = getLibraryNames(filepath)

	newid = None
	newdict = {}

	if os.path.basename(libnames[0]) in libfilenames:
		newid = os.path.join(libpath, os.path.basename(libnames[0]))
		if loader_path:
			newid = os.path.basename(libnames[0])

	for libname in libnames[1:]:
		if os.path.basename(libname) in libfilenames:
			newdict[libname] = os.path.join(libpath, os.path.basename(libname))
			if loader_path:
				relative_path = relativePath(os.path.dirname(filepath), libpath)
				newdict[libname] = os.path.join("@loader_path", relative_path, os.path.basename(libname))

	setLibraryNames(filepath, newdict, newid)

def setLibraryPathRecursive(dirpath, filepattern, libfilenames, libpath = "", loader_path = False):
	"""recursively update a set of library names be relative to a given
	directory"""

	if type(filepattern) != type(re.compile("")):
		filepattern = re.compile(filepattern)

	def dirfunc(arg, dirname, fnames):
		for fname in fnames:
			if filepattern.search(fname):
				setLibraryPath(os.path.join(dirname, fname), libfilenames, libpath, loader_path)

	os.path.walk(dirpath, dirfunc, None)

def relativePath(fromdir, todir):
	"""determine the relative path from one directory to another"""

	fromdir = os.path.abspath(fromdir)
	todir = os.path.abspath(todir)

	commonpath = os.path.commonprefix([fromdir, todir])

	fromdiff = fromdir[len(commonpath):]
	todiff = todir[len(commonpath):]

	return os.path.join(fromdiff.count(os.path.sep)*"../", todiff)

if __name__ == "__main__":

	setLibraryPathRecursive("rosetta", "\\.(so|dylib)$", ROSETTA_LIBS, os.path.abspath("rosetta"), True)
