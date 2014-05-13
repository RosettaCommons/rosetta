# /usr/bin/env python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os
import re
import sys
#import uuid
import random

proj = sys.argv[1].rstrip()
prefixdir = '../../src'
vsprefixdir = '../../../src'
sourcesFileName = '../../src/' + proj + '.src.settings'
vsFileName = '../VisualStudio/' + proj + '/' + proj + '.vcxproj'
vsFilterFileName = '../VisualStudio/' + proj + '/' + proj + '.vcxproj.filters'

vsLines = open(vsFileName, 'r').readlines()
vsPreLines, vsPostLines = [], []
addTo = 0
for vsLine in vsLines:
	if vsLine.find('<ItemGroup>') >= 0:
		addTo = 1
	elif vsLine.find('</ItemGroup>') >= 0 and addTo > 0:
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
vsFileLines += ['  <ItemGroup>\r\n']

vsIncludeLines = []
vsIncludeLines += ['  <ItemGroup>\r\n']

paths = {}

execfile(sourcesFileName)
for dir, srcfiles in sources.iteritems():
	if len(dir) != 0 and dir[-1] != '/':
		dir += '/'
	
	if not os.path.isdir(prefixdir + '/' + dir):
		print "Skipping non-existent directory:", dir
		continue	

	hdrfiles = os.listdir(prefixdir + '/' + dir)
	hdrfiles = [hdrfile for hdrfile in hdrfiles if hdrfile.find('.h') >= 0]

	srcfiles = [srcfile + '.cc' for srcfile in srcfiles if not srcfile.endswith('.cu')]

	allfiles = sorted(hdrfiles + srcfiles)

	for allfile in allfiles:
		allfile_path = dir + allfile
		allfile = vsprefixdir + '/' + dir + allfile
		allfile = allfile.replace('/', '\\')

		path = allfile_path[0:allfile_path.rfind('/')]
		if not paths.has_key(path):
			paths[path] = []
		paths[path] += [ allfile ]

		if allfile.find('.cc') == len(allfile) - 3:
			objfile = proj + '\\' + allfile[len(vsprefixdir) + 1:-3]
			objfile = objfile.replace('/', '__')
			objfile = objfile.replace('\\', '__')

			vsFileLines += ['    <ClCompile Include="'+allfile+'">\r\n']
			for config in ['BoincDebug|Win32','BoincRelease|Win32','Debug|Win32','Release|Win32']:
				vsFileLines += ['      <ObjectFileName Condition="\'$(Configuration)|$(Platform)\'==\'' + config + '\'">$(IntDir)' + objfile + '.obj</ObjectFileName>\r\n']
				vsFileLines += ['      <XMLDocumentationFileName Condition="\'$(Configuration)|$(Platform)\'==\'' + config + '\'">$(IntDir)' + objfile + '.xdc</XMLDocumentationFileName>\r\n']
			vsFileLines += ['    </ClCompile>\r\n']
			
		else:
			vsIncludeLines += ['    <ClInclude Include="' + allfile + '" />\r\n']

vsFileLines += ['  </ItemGroup>\r\n']
vsIncludeLines += ['  </ItemGroup>\r\n']

lines = vsPreLines + vsFileLines + vsIncludeLines + vsPostLines
open(vsFileName, 'w').writelines(lines)

# Filters
vsFilterLines = ['<?xml version="1.0" encoding="utf-8"?>\r\n']
vsFilterLines += ['<Project ToolsVersion="4.0" xmlns="http://schemas.microsoft.com/developer/msbuild/2003">\r\n']
vsFilterLines += ['  <ItemGroup>\r\n']
for path in paths.keys():
	filter_tag = '%255c' + path.replace('/', '%255c') + '%255c'
	#id = str(uuid.uuid5(uuid.uuid1(), path))
	id = "%06X-%04X-%04X-%04X-%012X" % (random.randint(0,0xFFFFFF), random.randint(0,0xFFFF), random.randint(0,0xFFFF), random.randint(0,0xFFFF), random.randint(0,0xFFFFFFFFFFFF))
	vsFilterLines += ['    <Filter Include="' + filter_tag + '">\r\n']
	vsFilterLines += ['      <UniqueIdentifier>{' + id + '}</UniqueIdentifier>\r\n']
	vsFilterLines += ['    </Filter>\r\n']
vsFilterLines += ['  </ItemGroup>\r\n']

vsFilterLines += ['  <ItemGroup>\r\n']
for path in paths.keys():
	for obj in paths[path]:
		filter_tag = '%255c' + path.replace('/', '%255c') + '%255c'
		vsFilterLines += ['    <ClCompile Include="' + obj + '">\r\n']
		vsFilterLines += ['      <Filter>' + filter_tag + '</Filter>\r\n']
		vsFilterLines += ['    </ClCompile>\r\n']
vsFilterLines += ['  </ItemGroup>\r\n']
vsFilterLines += ['</Project>\r\n']
open(vsFilterFileName, 'w').writelines(vsFilterLines)
