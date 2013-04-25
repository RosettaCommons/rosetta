#! /usr/bin/env python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os
import re
import sys

proj = sys.argv[1].rstrip()
prefixdir = '../../src'
vsprefixdir = '../../../src'
sourcesFileName = '../../src/' + proj + '.src.settings'
vsFileName = '../VisualStudio/' + proj + '/' + proj + '.vcproj'

vsLines = open(vsFileName, 'r').readlines()
vsPreLines, vsPostLines = [], []
addTo = 0
for vsLine in vsLines:
	if vsLine.find('<Files>') >= 0:
		addTo = 1
	elif vsLine.find('</Files>') >= 0:
		addTo = 2
	else:
		if addTo == 0:
			vsPreLines.append(vsLine)
		if addTo == 2:
			vsPostLines.append(vsLine)

# Additional include directories
pattern = re.compile('AdditionalIncludeDirectories="(.*)"')
for (i, line) in enumerate(vsPreLines):
	match = pattern.search(line)
	if match:
		setting = match.group(1) + ';..\..\..\external\dbio;'
		vsPreLines[i] = 'AdditionalIncludeDirectories="%s"\n' % setting

vsFileLines = []
vsFileLines += ['\t<Files>\r\n']

execfile(sourcesFileName)
for dir, srcfiles in sources.iteritems():
	if len(dir) != 0 and dir[-1] != '/':
		dir += '/'

	vsFileLines += ['\t\t<Filter\r\n']
	vsFileLines += ['\t\t\tName="\\' + dir.replace('/', '\\') + '"\r\n']
	vsFileLines += ['\t\t\t>\r\n']
	
	hdrfiles = os.listdir(prefixdir + '/' + dir)
	hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.find('.h') >= 0]

	srcfiles = [srcfile + '.cc' for srcfile in srcfiles if not srcfile.endswith('.cu')]

	allfiles = sorted(hdrfiles + srcfiles)

	allfiles = [vsprefixdir + '/' + dir + allfile for allfile in allfiles]
	allfiles = [allfile.replace('/', '\\') for allfile in allfiles]

	for allfile in allfiles:
		vsFileLines += ['\t\t\t<File\r\n']
		vsFileLines += ['\t\t\t\tRelativePath="' + allfile + '"\r\n']
		vsFileLines += ['\t\t\t\t>\r\n']

		if allfile.find('.cc') == len(allfile) - 3:
			objfile = proj + '\\' + allfile[len(vsprefixdir) + 1:-3]
			objfile = objfile.replace('/', '__')
			objfile = objfile.replace('\\', '__')

			for cfgName in ['Debug|Win32', 'Release|Win32', 'BoincDebug|Win32' , 'BoincRelease|Win32']:
				vsFileLines += ['\t\t\t\t<FileConfiguration\r\n']
				vsFileLines += ['\t\t\t\t\tName="' + cfgName + '"\r\n']
				vsFileLines += ['\t\t\t\t\t>\r\n']
				vsFileLines += ['\t\t\t\t\t<Tool\r\n']
				vsFileLines += ['\t\t\t\t\t\tName="VCCLCompilerTool"\r\n']
				vsFileLines += ['\t\t\t\t\t\tObjectFile="$(IntDir)\\' + objfile + '.obj"\r\n']
				vsFileLines += ['\t\t\t\t\t\tXMLDocumentationFileName="$(IntDir)\\' + objfile + '.xdc"\r\n']
				vsFileLines += ['\t\t\t\t\t/>\r\n']
				vsFileLines += ['\t\t\t\t</FileConfiguration>\r\n']

		vsFileLines += ['\t\t\t</File>\r\n']

	vsFileLines += ['\t\t</Filter>\r\n']

vsFileLines += ['\t</Files>\r\n']

lines = vsPreLines + vsFileLines + vsPostLines
open(vsFileName, 'w').writelines(lines)
