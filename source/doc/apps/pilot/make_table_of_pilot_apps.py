#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--file", dest="filename", default="../../../src/pilot_apps.src.settings.all",
  action="store", type="string", metavar="FILE",
  help="Name of pilot_apps.src.settings.all file from which to generate table of applications and demos" )

parser.add_option("-d", "--doxygen_file", dest="doxyfile", default="PilotAppsInfo.dox",
  action="store", type="string", help="Name of Dxygen file to write code to", metavar="STRING")

(options, args) = parser.parse_args()



f = open( options.filename, "r")
d = open( options.doxyfile, "w")


d.write("/// @page pilot_apps_info Pilot Application Information Page\n")
d.write("///\n")
d.write("///This page is automatically generated using the make_table_of_pilot_apps.py.\n")
d.write("///The script reads in the pilot_apps.src.settings file. It looks for line with \/\* \*\/ comments.\n")
d.write("///It splits the string inside by commas. The first three fields are included in the table\n")
d.write("///\n")
d.write("///All pilot applications should have entries in this page.\n")
d.write("///If you add a new pilot appication please add the information and regenerate this doxygen page using the above mentioned script.\n")
d.write("///Call python make_table_of_pilot_apps.py -h for options\n")
d.write("///\n")
d.write("/// <table border=\"1\">\n")
d.write("/// <tr>\n")
d.write("/// <td> Application </td>\n")
d.write("/// <td> Developers </td>\n")
d.write("/// <td> Documentation Location</td>\n")
d.write("/// <td> Demo File Location</td>\n")
d.write("/// </tr>\n")

for line in f :
  begin = line.find('/*')
  end = line.find('*/')
  start_name = line.find('"')
  end_name = line.find('"', start_name +1 )
  if begin != -1 and end != -1 :
    line_info = line[begin+1:end-1]
    d.write("/// <tr>\n")
    name = line[start_name+1:end_name]
    d.write("/// <td> " + name + " </td>\n")
    app_info = line_info.split(",")
    for i in range(0,3) :
      d.write("/// <td>" + app_info[i] + "</td>\n")
    d.write("/// </tr>\n")

d.write("/// </table>\n")
