#Author: 
 Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#About: 
 This code generates pre-configured files for Rosetta classes and files.  It can be run from any directory.
 Any paths needed will be created, templates are placed within source. User name and email are parsed from git.

#Running the Code

##Running the src generator:

 Use ./generate_templates.py --help for a full listing.

Example: 
 
    ./generate_templates.py --type general_class --class_name TestTemplates --brief "a class for testing templates" --namespace protocols testing


 This will generate hh, cc, fwd.hh, and Creator templates within src/protocols/testing configured for the TestTemplates name and inheriting from ref count.


##Running the Unit Test generator:
 Use ./generate_unit_test_template.py --help for a full listing.

 Example: 
    
    ./generate_unit_test_templates.py --class_name TemplateTests --brief "a unit test for testing templates" --outdirs protocols templates
  
 This will generate the unit test hh file in src/test/protocols/templates with the UnitTest class being TemplateTests.


##Running the Application generator:
 Use ./generate_app_template_JD2.py --help for a full listing.

 Example without opts: 
 
    ./generate_app_template_JD2.py --pilot --user_name jadolfbr --app_name test_templates --class_name TemplateTest --brief "an app for testing templates"


 _A more complicated example is giving a list of relevant options for JD2 through --app_options and optionally specificying new LOCAL options._
 _in::file::s and in::file::l are automatically added to relevant options. Useful for pilot apps_


 Example with relavent and local opts: 
 
     ./generate_app_template_JD2.py --pilot --user_name jadolfbr --app_name test_templates
      --class_name TemplateTest --brief "An app for testing templates" --app_options group1::bool1 group2::bool2 group1::group2::int1
      --boolean_opt group1::bool1 group2::bool2 --integer_opt group1::group2::int1



#Dev Information
##Adding Templates
Adding new src code templates:
 Make a new directory within the src directory and populate with files.  Names will be replaced.
 If Creator is in the name of the template, will concatonate that at the end.
 This will become a new type for one to pass as an option.  See other templates for examples.  Use replacements given below.
 To add a new replacement type, add to python code in self.replacements that returns a string and add a description to this file.


##Current Replacement Variables:

###General:
 --name-- : User name
 
 --email-- : User email
 
 --class-- : Class name

 --brief-- : The class or file brief


###Path Matching:
--path-- : Relative path to the file.  Will use namespace or dir_override option cmd-line option to create string. (ex: protocols/analysis)
 
--path_underscore-- :relative path to file for ifdefs

###Namespace Matching:
--namespace-- : The encapsulating namespace block of code.

--namespace_dot-- : The namespace with dots for Tracer objects

--end_namespace-- : The ending namespace block to finish encapsulation.

###Template Type Specific Matches

####Residue Selector
--res_sel_creator-- : Places the residue selector creator path in the template if namespace is not core
For ResidueSelectors, Creators are in core are contained in one one file.

####Unit Tests
--test_functions-- : Setup the names of the test unit test function if option is given.

####JD2 Application
--app_name-- : Name of the application file

--app_options-- : Block of options.register options.

--new_app_options_out-- : Block outside of try for new local options - OPT_KEY(type, name) or OPT_1GRP_KEY(type, namespace, name) etc.

--new_app_options_in-- : Block inside of try for new local options  - NEW_OPT(namespace_and_name, brief, default) etc.

