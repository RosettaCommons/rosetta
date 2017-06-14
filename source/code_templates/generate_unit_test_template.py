#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   generate_templates.py
## @brief  Script for generating Rosetta files for classes and files
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#See Readme for use.
from __future__ import print_function

import os
work_dir = os.getcwd()
if os.path.dirname(__file__): os.chdir(os.path.dirname(__file__))


from argparse import ArgumentParser

import glob
import sys
import re
from generate_templates import GenerateRosettaTemplates

class GenerateUnitTestTemplate(GenerateRosettaTemplates):
    """
    Template Generator specific for Unit Tests
    """
    def __init__(self):

        parser = ArgumentParser(description="This class is used to generate Rosetta templates for use in any IDE. "
                                            "See the README for more detailed instructions.  ")

        required = parser.add_argument_group("Required")

        required.add_argument("--class_name", "-c",
                            help = "The name of the test class.  Should end with Test or Tests")

        required.add_argument("--brief", "-b",
                            help = "A brief description of the test.  Enclose in quotes.")

        required.add_argument("--outdirs",
                            help = "Where to place test file. Separate by spaces. "
                                   "Not used IN file."
                                   "Ex: --out protocols antibody",
                            nargs='*',
                            default = [])


        optional = parser.add_argument_group("Optional")

        optional.add_argument("--test_functions",
                            help = "Optional list of test name functions.",
                            nargs='*',
                            default = [])



        GenerateRosettaTemplates.__init__(self, "unit_test", parser)
        self.options.namespace = self.options.outdirs

        self.replacement["--test_functions--"] = lambda: self.get_test_function_block()

    def get_base_outdir(self):
        if self.options.test:
            return os.path.join(work_dir,"test_unit_test")
        else:
            return self.test_dir

    def get_test_function_block(self):

        if not self.options.test_functions:
            return "\n"
        else:
            blocks = []
            for func_name in self.options.test_functions:
                if not re.search("test", func_name):
                    func_name = "test_"+func_name
                block = "\n\tvoid "+func_name+"(){\n\n\n\t}"
                blocks.append(block)
            return "\n".join(blocks)

    def print_dev_help(self):
        print("\ngit add "+os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())+"/"+\
                      self.get_option("class_name", fail_on_none=False)+".cxxtest.hh\n")


if __name__ == "__main__":
    generate_files = GenerateUnitTestTemplate()
    generate_files.apply()


