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
from collections import defaultdict
from generate_templates import GenerateRosettaTemplates

class GenerateAppTemplate(GenerateRosettaTemplates):
    """
    Template Generator specific for applications using JD2
    """
    def __init__(self):

        self.types = sorted(
            [os.path.basename(d) for d in glob.glob(os.path.join("application", "*")) if os.path.isdir(d)])

        option_types = ["Boolean", "Integer", "String", "Real"]
        option_types.extend([x+"Vector" for x in option_types])


        self.option_type_names = defaultdict() #Option types to the names of the actual options.


        parser = ArgumentParser(description="This class is used to generate Rosetta templates for use in any IDE. "
                                            "See the README for more detailed instructions.  ")

        required = parser.add_argument_group("Required")

        required.add_argument("--type",
                            help = "The type of template you will be needing.",
                            required = True,
                            choices = self.types)

        required.add_argument("--app_name", "-a",
                            help = "The name of the app.",
                            required = True)



        required.add_argument("--mover_name", "-c",
                            help = "The name of the main Mover you are calling (JD2 especially).",
                            required = True)

        required.add_argument("--brief", "-b",
                            help = "A brief description of the app.  Enclose in quotes.",
                            required = True)

        required.add_argument("--mover_namespace",
                            help = "Main Mover namespace to add. Will add this hh file for include",
                            nargs='*',
                            default = [])



        optional = parser.add_argument_group("Optional")

        optional.add_argument("--pilot", "-p",
                            help = "Signify that this is a pilot app",
                            default = False,
                            action = "store_true")

        optional.add_argument("--user_name",
                            help = "User name if Pilot app")

        optional.add_argument("--app_dir",
                            help = "Any app DIR if public app or directory in pilot app user name directory",
                            nargs = '*',
                            default = [] )


        optional.add_argument("--app_options",
                            help = "Register needed app options. "
                                   "in:file:s and in:file:l are used by default",
                            nargs="*",
                            default=["in::file::s", "in::file::l"])


        new_opts = parser.add_argument_group("Optional list of app options. (not options_rosetta.py opts)  Ex: antibody::graft_L1  "
                                             "(Only recommended for Pilot Apps)")


        for arg_type in option_types:
            #arg_type = arg_type.lower()
            self.option_type_names[arg_type] = arg_type.lower()+"_opt"

            new_opts.add_argument("--"+arg_type.lower()+"_opt",
                                  help = "Optional List of local new "+arg_type+" options with namespaces.  Ex: grid::grid_length",
                                  nargs = '*',
                                  default = [])


        GenerateRosettaTemplates.__init__(self, "application", parser)

        self.options.class_name = self.options.mover_name #Same thing, different name.

        ##Extend Matching
        self.replacement["--app_name--"] = lambda: self.get_option("app_name", fail_on_none=True)
        self.replacement["--app_options--"] = lambda: self.get_app_options()
        self.replacement["--new_app_options_out--"] = lambda: self.get_new_app_options_out()
        self.replacement["--new_app_options_in--"] = lambda: self.get_new_app_options_in()
        self.replacement["--mover_namespace--"] = lambda: self.get_mover_namespace()
        self.replacement["--mover_path--"] = lambda: self.get_mover_path()

    def get_mover_namespace(self):
        if self.options.mover_namespace:
            return "::".join(self.options.mover_namespace)
        else:
            return ""

    def get_mover_path(self):
        if not self.options.mover_namespace:
            return "\n"
        else:
            return "#include <"+"/".join(self.options.mover_namespace)+"/"+self.options.class_name+".hh>\n"

    def get_base_outdir(self):
        if self.options.test:
            return os.path.join(work_dir,"test_src")
        else:
            return self.src_dir

    def get_outname(self):
        return self.get_option("app_name", fail_on_none=True)


    def get_app_options(self):
        """
        Get a string for registering required app options.
        :rtype: str
        """
        default_opts = ["in::file::s", "in::file::l"]
        for opt in default_opts:
            if not opt in self.options.app_options:
                self.options.app_options.append(opt)

        lines = ["\n"]
        line_fmt = "\toption.add_relevant( {opt} );"
        for opt in self.options.app_options:

            line = line_fmt.format(opt=opt)
            lines.append(line)
        return "\n".join(lines)

    def get_new_app_options_out(self):
        """
        Get a string for registering new app options if not empty (Outside Try).
        :rtype: str
        """
        opt_lines = ["\n"]
        line_fmt = "OPT_{opt_n}KEY( {opt_type}, {opts} )"
        all_opts = vars(self.options)
        for opt_type in self.option_type_names:
            opt_name = self.option_type_names[opt_type]
            if all_opts.has_key(opt_name) and all_opts[opt_name]:
                for opt in all_opts[opt_name]:
                    optSP = opt.split("::")
                    if len(optSP) == 1:
                        line = line_fmt.format(opt_type=opt_type, opts=opt, opt_n="")
                    else:
                        n_group = len(optSP) -1
                        opt_n = str(n_group)+"GRP"
                        grouped = ", ".join(optSP)
                        line = line_fmt.format(opt_type=opt_type, opts=grouped, opt_n=opt_n)
                    opt_lines.append(line)

        return "\n".join(opt_lines)


    def get_new_app_options_in(self):
        """
        Get a string for setting new app options if not empty (Inside Try).
        :rtype str:
        """

        opt_lines = ["\n//this won't compile until you fill in brief and default yourself"]
        line_fmt = "\t\tNEW_OPT( {opt_name}, brief, default );"
        all_opts = vars(self.options)
        for opt_type in self.option_type_names:
            opt_name = self.option_type_names[opt_type]
            if all_opts.has_key(opt_name) and all_opts[opt_name]:
                for opt in all_opts[opt_name]:
                    line = line_fmt.format(opt_name=opt)
                    opt_lines.append(line)
        return "\n".join(opt_lines)



    def get_path_list(self):
        """
        Get a list of the path directory hierarchy.
        :rtype: list
        """

        if self.options.pilot:
            if not self.options.user_name:
                sys.exit("Please give the option --user_name for pilot apps")
            else:
                l = ["apps", "pilot", self.options.user_name]
                if self.options.app_dir:
                    l.extend(self.options.app_dir)

                return l
        else:
            l = ["apps", "public"]
            if self.options.app_dir:
                l.extend(self.options.app_dir)

            return l


if __name__ == "__main__":
    generate_files = GenerateAppTemplate()
    generate_files.apply()


