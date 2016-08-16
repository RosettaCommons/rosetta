#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/rosetta_tools/SetProtocolOptions3-3.py
## @brief  Used originally for manually currated options and descriptions.  
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


import pickle
import os
import re
#This is the file used to generate the protocol options binary.
#Use - Python SetProtocolOptions3-3.py
OPTIONS = dict(); #List of options for each application and their descriptions
APPLICATIONS = dict(); #List of Applications and a description for each known one.

def location():
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported from anywhere.
        """
        
        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        print pathSP
        return pathSP
def readInOptionData(pathToData, DICTIONARY):
      directorylist =os.listdir(pathToData)
      for f in directorylist:
            if re.search(".txt", f) and (re.search("\._", f)== None) and (re.search("~", f)== None):#This shows how stupid python/linux can be...
                  print "Reading Application in: "+f
                  FILE = open(pathToData+"/"+f, 'r')
                  fileSP = f.split("."); app = fileSP[0]
                  DICTIONARY[app]=dict()
                  for line in FILE:
			if line == "\n":
				break
                        lineSP = line.split("==")
                        DICTIONARY[app][lineSP[0]]=lineSP[1]
                  FILE.close()
      return DICTIONARY
                                 



#Descriptions of each of the Applications for which there is sufficient documentation.
#Most descriptions come from rosettacommons.org






if __name__ == '__main__':
            DescriptionDIR = "/AppDescriptions/Rosetta3-3"
	    OptionDIR = "/AppOptions/Rosetta3-3"
            pwd = location()
	    print "Reading Description Data..."
            APPLICATIONS = readInOptionData(pwd[0]+DescriptionDIR, APPLICATIONS)
            #print APPLICATIONS
	    print "Reading Option Data..."
	    OPTIONS= readInOptionData(pwd[0]+OptionDIR, OPTIONS)
            #print OPTIONS
            pickle.dump(OPTIONS, open("Rosetta3-3.p", 'wb'))
            pickle.dump(APPLICATIONS, open("Rosetta3-3Apps.p", 'wb'))
