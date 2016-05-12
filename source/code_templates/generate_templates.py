#!/usr/bin/env python
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   generate_templates.py
## @brief  Script for generating Rosetta files for classes and files
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#See Readme for use.

import os
work_dir = os.getcwd()

if os.path.dirname(__file__): os.chdir(os.path.dirname(__file__))
from argparse import ArgumentParser

import glob
import sys
import re
import subprocess

template_types = ["src", "application", "unit_test"]
residue_selector_namespace = ["core", "select", "residue_selector"]


class GenerateRosettaTemplates(object):
    def __init__(self, template_type, parser):

        testing_args = parser.add_argument_group("Testing")
        testing_args.add_argument("--test",
                            help = "Indicate that we are running in test mode.  "
                                   "All output files will go to the current work dir.",
                            default = False,
                            action = "store_true",
                            required = False)

        if len(sys.argv) < 2:
            #parser.print_usage()
            parser.print_help()
            sys.exit()


        self.options = parser.parse_args()
        self.template_type = template_type
        self.start_pwd = os.getcwd()

        #os.chdir(template_type)

        self.source_dir = os.path.abspath("..")
        self.src_dir = os.path.abspath("../src")
        self.test_dir = os.path.abspath("../test")


        #Functions to run for replacement.  Must return a string.

        self.replacement = {
            "--name--" : lambda: self.get_name(),
            "--email--": lambda: self.get_email(),
            "--class--": lambda: self.get_option("class_name", fail_on_none=True),
            "--brief--": lambda: self.get_option("brief", fail_on_none=True),
            "--path--": lambda: self.get_outfile_rel_path(),
            "--path_underscore--": lambda: self.get_path_underscore(),
            "--namespace--": lambda: self.get_namespace_lines(),
            "--namespace_dot--": lambda: self.get_namespace_char_concat("."),
            "--namespace_colon--": lambda: self.get_namespace_char_concat(":"),
            "--end_namespace--": lambda: self.get_end_namespace(),
            "--res_sel_creator--": lambda: self.get_res_sel_creator_path()
        }

        #### General Template Paths ###
        self.fwd_class_template = os.path.join(template_types[0], "RosettaClass.fwd.hh")

    def apply(self):

        if hasattr(self.options, "type") and self.options.type:
            template_dir = os.path.join(self.template_type, self.options.type)
            files = glob.glob(os.path.join(template_dir, '*'))
        else:
            template_dir = self.template_type
            files = glob.glob(os.path.join(template_dir, "*"))

        print "Type: "+self.template_type

        #Add forward class template for most types
        if self.fwd_class_template not in files \
                and self.template_type == template_types[0]\
                and self.options.type != "util":
            files.append(self.fwd_class_template)

        #Special case residue selector creator
        if hasattr(self.options, "type") and self.options.type == "residue_selector" and \
            self.options.namespace == residue_selector_namespace:

            files = [f for f in files if os.path.basename(f) != "ResidueSelectorCreator.hh"]



        outdir = os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        matches = self.replacement.keys()

        for template in sorted(files):
            extension = "."+".".join(os.path.basename(template).split(".")[1:])

            if hasattr(self.options, "type") and self.options.type == "util":
                out_basename = "util" + extension
            elif re.search("Creator", os.path.basename(template)):
                out_basename = self.get_outname()+"Creator"+extension
            else:
                out_basename = self.get_outname()+extension

            new_file_name = os.path.join(outdir, out_basename)

            print "Creating "+new_file_name

            INFILE = open(os.path.abspath( template ), 'r')
            out_lines = []
            for line in INFILE:
                newline = line
                for key in matches:
                    if re.search(key, line):
                        newline = newline.replace(key, self.replacement[key]())

                out_lines.append(newline)

            INFILE.close()


            #Write the output file ONLY if all things are OK.
            OUTFILE = open(new_file_name, 'w')
            for line in out_lines:
                OUTFILE.write(line)

            OUTFILE.close()


        self.print_dev_help()


    def get_base_outdir(self):
        if self.options.test:
            return os.path.join(work_dir,"test_src")
        else:
            return self.src_dir

    def print_dev_help(self):
        """
        Printed at the end of template generation to help devs know where things are supposed to go, etc.
        """



        if hasattr(self.options, "type"):

            if self.options.type == "util":
                print "\ngit add "+os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())+"/util"+"*"
            else:
                print "\ngit add "+os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())+"/"+\
                      self.get_option("class_name", fail_on_none=False)+"*"


            if self.options.type == "residue_selector" and self.options.namespace == residue_selector_namespace:
                print "\nA Creator class should be declared in core::select::residue_selector::ResidueSelectorCreators.hh"

            if self.options.type == "mover":
                print "\nMover Creator should be registered in (protocols.7) \n" \
                      "   protocols/init/init.MoverRegistrators.ihh and \n" \
                      "   protocols/init/init.MoverCreators.ihh\n"



    ######## Replacement Functions#############
    def get_option(self, option_name, fail_on_none = True):
        opts = vars(self.options)
        if not opts.has_key(option_name) and fail_on_none:
            sys.exit(option_name+" is necessary.  Pass it as an argument.")
        elif not opts[option_name] and fail_on_none:
            sys.exit(option_name+" is necessary.  Pass it as an argument.")
        elif opts.has_key(option_name):
            return opts[option_name]
        else:
            return None

    def get_name(self):
        return subprocess.check_output("git config user.name", shell=True).strip()

    def get_email(self):
        return subprocess.check_output("git config user.email", shell=True).strip()

    def get_outname(self):

        return self.get_option("class_name", fail_on_none= True)


    def get_outfile_rel_path(self):
        """
        Get the rel path line.  Ex: protocols/antibody
        :rtype: str
        """
        return "/".join(self.get_path_list())

    def get_path_underscore(self):
        """
        Get the path with underscore for ifdefs
        :rtype: str
        """
        return "_".join(self.get_path_list())

    def get_namespace_lines(self):
        """
        Get the namespace declaration lines for hh and cc file
        :rtype: str
        """
        return "\n".join(["namespace "+n+" {" for n in self.get_option("namespace", fail_on_none=True) ] )


    def get_namespace_char_concat(self, char):
        return char.join(self.get_option("namespace", fail_on_none=True))

    def get_end_namespace(self):
        """
        Get the end of the namespace declaration at end of file.
        :rtype: str
        """
        return "\n".join(["} //"+n for n in self.get_option("namespace", fail_on_none=True)])

    def get_path_list(self):
        """
        Get a list of the path directory hierarchy.
        :rtype: list
        """
        if hasattr(self.options, "dir_override") and self.options.dir_override:
            return self.options.dir_override
        elif hasattr(self.options, "namespace") and self.options.namespace:
            return self.options.namespace
        else:
            sys.exit("Path not defined.  Either set the path override or pass a namespace")

    def get_res_sel_creator_path(self):
        """
        Places the residue selector creator path in the template if namespace is not core
        For ResidueSelectors, Creators are in core are contained in one one file.
        :rtype: str
        """
        if self.options.namespace == residue_selector_namespace:
            return "<core/select/residue_selector/ResidueSelectorCreators.hh>\n"
        else:
            return "<"+self.replacement["--path--"]()+"/"+self.replacement["--class--"]()+"Creator.hh>\n"



class GenerateGeneralTemplates(GenerateRosettaTemplates):
    """
    Template Generator specifically for general rosetta classes and files (movers, task_ops, etc.)
    """
    def __init__(self, template_type_name = template_types[0]):
        self.types = sorted( [os.path.basename(d) for d in glob.glob(os.path.join(template_type_name,"*")) if os.path.isdir(d)] )
        print "Found template types: "+repr(self.types)


        parser = ArgumentParser(description="This class is used to generate Rosetta templates for use in any IDE. "
                                            "See the README for more detailed instructions.  ")

        required = parser.add_argument_group("Required")
        required.add_argument("--type",
                            help = "The type of template you will be needing.",
                            required = True,
                            choices = self.types)

        required.add_argument("--class_name", "-c",
                            help = "The name of the class you are creating if not creating util files.")

        required.add_argument("--brief", "-b",
                            help = "A brief description of the class/file.  Enclose in quotes.",
                            required = True)

        required.add_argument("--namespace",
                            help = "Namespace needed for file.  Separate by spaces. "
                                   "Default is to place files in this directory. "
                                   "Ex: --namespace protocols antibody",
                            nargs='*',
                            default = [],
                            required = True)

        optional = parser.add_argument_group("Optional")

        optional.add_argument("--dir_override",
                            help = "List of dir names if path to output file is different than namespace.",
                            nargs='*',
                            default = [])

        GenerateRosettaTemplates.__init__(self, template_type_name, parser)




if __name__ == "__main__":

    generate_files = GenerateGeneralTemplates()
    generate_files.apply()


