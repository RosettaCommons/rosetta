#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rosie/__main__.py
## @brief  main file for ROSIE-2 scripts package
## @author Sergey Lyskov

import os, sys, importlib

from argparse import ArgumentParser


#
# list of ROSIE scipt-modules, add your app name here and put related code int <app-name>.py file
#
script_modules = '''
util docking stabilize-pm stabilize-cluster
'''

def rosie_scripts():
    ''' return map script-name --> excution function
    '''
    scripts = {}
    for module_name in script_modules.split():
        module = importlib.import_module('rosie.' + module_name)
        scripts.update( module.register() )

    return scripts




def execute_flag_file(file_name):
    ''' parse and "execute" given flag-file
    '''
    scripts = rosie_scripts()
    with open(file_name) as f:
        for line in f.read().split('\n'):
            options = line.split()
            if options:
                print(line)
                command = options.pop(0)

                if command in scripts: scripts[command](*options)
                else:
                    print(f'ERROR: unknown command {command!r}, terminating...')
                    sys.exit(1)




def main(args) -> None:
    ''' ROSIE scripts '''

    parser = ArgumentParser(description=main.__doc__)

    parser.add_argument('-w', '--working-directory', default=os.getcwd(), help="Specify scripts working directory")

    parser.add_argument("flags")

    options = parser.parse_args()

    print( f'working directory: {options.working_directory}, flags file: {options.flags}' )

    os.chdir(options.working_directory)

    execute_flag_file(options.flags)




if __name__ == "__main__": main(sys.argv)
