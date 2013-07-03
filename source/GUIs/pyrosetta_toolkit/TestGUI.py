#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/TestGUI.py
## @brief  Tests the GUI using Sergey's TestBindings.py 
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import sys
import os


def get_py_files(dir_):
    """
    Sergey Lyskov
    """
    return [dir_ + '/' + f for f in os.listdir(dir_) if f.endswith('.py')]



if __name__ == '__main__':
    """
    To run this test script: Supply full path to TestBindings.py.  It is in root PyRosetta directory / or in trunk: src/python/bindings
    """
    test_bindings = os.path.abspath("../../src/python/bindings")+"/TestBindings.py"
    
    #Locate TestBindings.
    if not os.path.exists(test_bindings):
        test_bindings = os.path.abspath("../../TestBindings.py")
        if not os.path.exists(test_bindings):
            print "Default location for test bindings not found.  Looking at script argument"
            test_bindings = sys.argv[0]
            if not os.path.exists(test_bindings):
                os.sys.exit("Could not locate TestBindings.py.  Please pass full path to test script.")
    
    tests = get_py_files("tests")
    os.system("python "+test_bindings+" "+" ".join(tests))