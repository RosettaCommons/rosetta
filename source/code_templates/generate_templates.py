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

from __future__ import absolute_import
import os
work_dir = os.getcwd()
# change the directory to this file's directory, if it is run from elsewhere.
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
            "--namespace_2colon--": lambda: self.get_namespace_char_concat("::"),
            "--namespace_underscore--": lambda: self.get_namespace_char_concat("_"),
            "--end_namespace--": lambda: self.get_end_namespace(),
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

        print("Type: ", self.template_type)

        #Add forward class template for most types
        if self.fwd_class_template not in files \
                and self.template_type == template_types[0]\
                and self.options.type != "util":
            files.append(self.fwd_class_template)

        #Special case residue selector creator
        #if hasattr(self.options, "type") and self.options.type == "residue_selector" and \
        #    self.options.namespace == residue_selector_namespace:
        #
        #    files = [f for f in files if os.path.basename(f) != "ResidueSelectorCreator.hh"]



        outdir = os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        matches = list(self.replacement.keys())

        for template in sorted(files):
            extension = "."+".".join(os.path.basename(template).split(".")[1:])

            if re.search("Creator", os.path.basename(template)):
                out_basename = self.get_outname()+"Creator"+extension
            else:
                out_basename = self.get_outname()+extension

            new_file_name = os.path.join(outdir, out_basename)

            print("Creating ",new_file_name)

            INFILE = open(os.path.abspath( template ), 'r')
            out_lines = []
            for line in INFILE:
                newline = line
                for key in matches:
                    if re.search(key, line):

                        newline = newline.replace(key, str(self.replacement[key]()))

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
            print("\nRemember to add your newly created files to the git repository:")
            print("\ngit add "+os.path.join(self.get_base_outdir(), self.get_outfile_rel_path())+"/"+\
                      self.get_option("class_name", fail_on_none=False)+"*")


            if self.options.type == "residue_selector":

                if self.options.namespace[0] == "core":
                    print("\nRegister in (core.6): \n"+\
                            "   "+self.get_base_outdir()+"/"+"core/init/init.cc")
                else:
                    print("\nRegister in (protocols.8):\n" \
                            "   "+self.get_base_outdir()+"/"+"protocols/init/init.ResidueSelectorCreators.ihh\n" \
                            "   "+self.get_base_outdir()+"/"+"protocols/init/init.ResidueSelectorRegistrators.ihh\n")

                print("\n See Wiki for How to serialize your Selector: \n" \
                    "   "+"https://wiki.rosettacommons.org/index.php/SerializationFAQ")

            elif self.options.type == "crosslinker_mover_helper" :
                print( "\n CrosslinkerMoverHelpers need to be registered in protocols_a.6 (where the CrosslinkerMoverHelper base class is registered) or higher." )
                print( " They also need to be tied into the CrosslinkerMover, since there is no factory system for helpers.  The src/protocols/cyclic_peptide/CrosslinkerMover.cc file must include the new header file, and CrosslinkerMover::crosslinkermover_helper_from_type() function's switch statement must be modified to create an instance of the new helper class.  It is also necessary to modify CrosslinkerMover::get_crosslinker_name(), CrosslinkerMover::provide_xml_schema() (particularly the linker_possibles string), and the Crosslinker enum class in the header file, CrosslinkerMover.hh.")

            elif self.options.type == "singleton":
                print("\n Singletons do not need to be registered.  " \
                      "To get the current instance from anywhere, use:\n" \
                      "  "+self.replacement["--namespace_2colon--"]()+"::"+self.replacement["--class--"]()+"::get_instance()->")

            elif self.options.type == "mover":
                print("\nMover Creator should be registered in (protocols.8) \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.MoverRegistrators.ihh and \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.MoverCreators.ihh\n")

            elif self.options.type == "features_reporter":
                print("\nFeature Reporter Creator should be registered in (protocols.8) \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.FeaturesReporterRegistrators.ihh and \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.FeaturesReporterCreators.ihh\n")
            elif self.options.type == "constraint_generator":
                print("\nConstraint Generator Creator should be registered in (protocols.8) \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.ConstraintGeneratorRegistrators.ihh and \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.ConstraintGeneratorCreators.ihh\n")

            elif self.options.type == "jd3_standard":
                     print("\n   "+"This template is for a standard JD3 app, however, much more complex apps can be created.  See the docs for more.\n")

            elif re.search("metric_", self.options.type):
                print("\nSimple Metric Creator should be registered in (protocols.8) \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.SimpleMetricRegistrators.ihh and \n" \
                      "   "+self.get_base_outdir()+"/"+"protocols/init/init.SimpleMetricCreators.ihh\n")

            elif re.search("jd3", self.options.type):
                print("\nDocs can be found out https://wiki.rosettacommons.org/index.php/JD3FAQ\n")

                if self.options.type == "jd3_standard_job_queen":
                    print("A tutorial can be found here: \n"
                          "https://www.rosettacommons.org/docs/wiki/development_documentation/tutorials/jd3_derived_jq/jd3_derived_jq_home\n")

                elif self.options.type == "jd3_job":
                    print("A tutorial can be found here: \n"
                          "https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/tutorial_job\n")

                elif self.options.type == "jd3_job_summary":
                    print("A tutorial can be found here: \n "
                          "https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/completed_job_summary\n")

    ######## Replacement Functions#############
    def get_option(self, option_name, fail_on_none = True):
        opts = vars(self.options)
        if option_name not in opts and fail_on_none:
            sys.exit(option_name+" is necessary.  Pass it as an argument.")
        elif not opts[option_name] and fail_on_none:
            sys.exit(option_name+" is necessary.  Pass it as an argument.")
        elif option_name in opts:
            return opts[option_name]
        else:
            return None

    def get_name(self):
        name = subprocess.check_output("git config user.name", shell=True).strip()
        if sys.version_info[0] > 2:
           name = name.decode(encoding="utf-8", errors="replace")
        return name

    def get_email(self):
        email = subprocess.check_output("git config user.email", shell=True).strip()
        if sys.version_info[0] > 2:
            email = email.decode(encoding="utf-8", errors="replace")
        return email

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
        
        return "\n".join(["} //"+n for n in self.get_option("namespace", fail_on_none=True)][::-1])

    def get_path_list(self):
        """
        Get a list of the path directory hierarchy.
        :rtype: list
        """
        if hasattr(self.options, "dir_override") and self.options.dir_override:

            #Catching things like protocols/antibody
            if len(self.options.dir_override) == 1:
                self.options.dir_override = self.options.dir_override[0]

                #Check to make sure someone did not try to use src/test to set directory.  Helps new users.
                if re.search('src', "/".join(self.options.dir_override)):
                    sys.exit("\nERROR. Please check your path.  Path is relative, src directory not needed. ( Ex - protocols/antibody )")

                #Not split.  We don't fail here, but instead give it as a list.
                if (re.search('/', self.options.dir_override)):
                    return self.options.dir_override.strip().split('/')
            else:
                return self.options.dir_override

        elif hasattr(self.options, "namespace") and self.options.namespace:
            return self.options.namespace
        else:
            sys.exit("Path not defined.  Either set the path/outdirs or pass a namespace")

    #def get_res_sel_creator_path(self):
    #    """
    #    Places the residue selector creator path in the template if namespace is not core
    #    For ResidueSelectors, Creators are in core are contained in one one file.
    #    :rtype: str
    #    """
    #    if self.options.namespace == residue_selector_namespace:
    #        return "<core/select/residue_selector/ResidueSelectorCreators.hh>\n"
    #    else:
    #        return "<"+self.replacement["--path--"]()+"/"+self.replacement["--class--"]()+"Creator.hh>\n"



class GenerateGeneralTemplates(GenerateRosettaTemplates):
    """
    Template Generator specifically for general rosetta classes and files (movers, task_ops, etc.)
    """
    def __init__(self, template_type_name = template_types[0]):
        self.types = sorted( [os.path.basename(d) for d in glob.glob(os.path.join(template_type_name,"*")) if os.path.isdir(d)] )
        print("Found template types: "+repr(self.types))


        parser = ArgumentParser(description="This class is used to generate Rosetta templates for use in any IDE. "
                                            "See the README for more detailed instructions.  ")

        required = parser.add_argument_group("Required")
        required.add_argument("--type",
                            help = "The type of template you will be needing: "+", ".join(self.types),
                            required = True,
                            metavar = "class_type",
                            choices = self.types)

        required.add_argument("--class_name", "-c",
                            help = "The name of the class you are creating. In case of util type, this will be the prefix for util files.")

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


