// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/rosetta_scripts/rosetta_scripts.cc
/// @brief The application file for rosetta_scripts, aka jd2_scripting or the parser
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added template script generation and "-parser:info" flag (in-app help for RosettaScripts).

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/BOINCJobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <core/types.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>

// Utility Headers

// Unit Headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>

#include <basic/Tracer.hh>

#include <fstream>

// Tracer
static THREAD_LOCAL basic::Tracer TR( "apps.public.rosetta_scripts.rosetta_scripts" );

// FUNCTION PROTOTYPES
void* my_main( void *);

/// @brief Prints out an empty template RosettaScript to the tracer.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void print_template_script();

/// @brief Prints out XSD information about the XML-accessible options for a given RosettaScipts-accessible
/// mover, filter, task operation, or residue selector.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool print_information( std::string const &component_name, std::stringstream &outstream );

/// @brief Prints out XSD information about the XML-accessible options for a given set of RosettaScipts-accessible
/// movers, filters, task operations, or residue selectors.
/// @details Calls the single string version.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void print_information( utility::vector1 < std::string > const &component_names );

/// @brief Saves the XSD to the given file.
void save_schema(  std::string const & filename );

// FUNCTION DECLARATIONS
void*
my_main( void *)
{
	protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it. That means that this can only be used with parser=true
	protocols::jd2::JobDistributor::get_instance()->go(mover);
	return 0 ;
}

/// @brief Prints out an empty template RosettaScript to the tracer.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
print_template_script() {
	TR << "No XML file was specified with the \"-parser:protocol <filename>\" commandline option.  In order for RosettaScripts to do something, it must be provided with a script." << std::endl;
	TR << "The following is an empty (template) RosettaScripts XML file:\n"
		<< "\n"
		<< "<ROSETTASCRIPTS>\n"
		<< "\t<SCOREFXNS>\n"
		<< "\t</SCOREFXNS>\n"
		<< "\t<RESIDUE_SELECTORS>\n"
		<< "\t</RESIDUE_SELECTORS>\n"
		<< "\t<TASKOPERATIONS>\n"
		<< "\t</TASKOPERATIONS>\n"
		<< "\t<FILTERS>\n"
		<< "\t</FILTERS>\n"
		<< "\t<MOVERS>\n"
		<< "\t</MOVERS>\n"
		<< "\t<APPLY_TO_POSE>\n"
		<< "\t</APPLY_TO_POSE>\n"
		<< "\t<PROTOCOLS>\n"
		<< "\t</PROTOCOLS>\n"
		<< "\t<OUTPUT />\n"
		<< "</ROSETTASCRIPTS>\n\n";
	TR << "At any point in a script, you can include text from another file using <xi:include href=\"filename.xml\" />." << std::endl;
	TR << "Variable substituion is possible from the commandline using the -\"parser:script_vars varname=value\" flag.  Any string of the pattern \"%%varname%%\" will be replaced with \"value\" in the script." << std::endl;
	TR << std::endl;
	TR << "The rosetta_scripts application will now exit." << std::endl;
	TR.flush();
}

/// @brief Prints out XSD information about the XML-accessible options for a given RosettaScipts-accessible
/// mover, filter, task operation, or residue selector.
/// @details Returns true for FAILURE to find the given component, false otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
print_information(
	std::string const &component_name,
	std::stringstream &outstream
) {
	bool missing( true );

	// 1. Check filters:
	protocols::filters::FilterFactory* filterfactory( protocols::filters::FilterFactory::get_instance() ); //Must be raw pointer; owning pointer does weird things with static singleton.
	if ( filterfactory->filter_exists( component_name ) ) {
		missing = false;
		outstream << "INFORMATION ABOUT FILTER \"" << component_name << "\":\n\n";
		utility::tag::XMLSchemaDefinition xsd;
		filterfactory->provide_xml_schema( component_name, xsd );
		outstream << xsd.human_readable_summary( component_name, "filter" );
	}

	// 2. Check movers:
	protocols::moves::MoverFactory* moverfactory( protocols::moves::MoverFactory::get_instance() ); //Must be raw pointer; owning pointer does weird things with static singleton.
	if ( moverfactory->mover_exists( component_name ) ) {
		if ( !missing ) outstream << "\n";
		missing = false;
		outstream << "INFORMATION ABOUT MOVER \"" << component_name << "\":\n\n";
		utility::tag::XMLSchemaDefinition xsd;
		moverfactory->provide_xml_schema( component_name, xsd );
		outstream << xsd.human_readable_summary( component_name, "mover" );
	}

	// 3. Check residue selectors:
	core::select::residue_selector::ResidueSelectorFactory* residueselectorfactory( core::select::residue_selector::ResidueSelectorFactory::get_instance() ); //Must be raw pointer; owning pointer does weird things with static singleton.
	if ( residueselectorfactory->has_type( component_name ) ) {
		if ( !missing ) outstream << "\n";
		missing = false;
		outstream << "INFORMATION ABOUT RESIDUE SELECTOR \"" << component_name << "\":\n\n";
		utility::tag::XMLSchemaDefinition xsd;
		residueselectorfactory->provide_xml_schema( component_name, xsd );
		outstream << xsd.human_readable_summary( component_name, "rs" );
	}

	// 4. Check task operations:
	core::pack::task::operation::TaskOperationFactory* taskfactory( core::pack::task::operation::TaskOperationFactory::get_instance() ); //Must be raw pointer; owning pointer does weird things with static singleton.
	if ( taskfactory->has_type( component_name ) ) {
		if ( !missing ) outstream << "\n";
		missing = false;
		outstream << "INFORMATION ABOUT TASK OPERATION \"" << component_name << "\":\n\n";
		utility::tag::XMLSchemaDefinition xsd;
		taskfactory->provide_xml_schema( component_name, xsd );
		outstream << xsd.human_readable_summary( component_name, "to" );
	}

	// 5. Check residue-level task operations:
	core::pack::task::operation::ResLvlTaskOperationFactory* rltaskfactory( core::pack::task::operation::ResLvlTaskOperationFactory::get_instance() ); //Must be raw pointer; owning pointer does weird things with static singleton.
	if ( rltaskfactory->has_type( component_name ) ) {
		if ( !missing ) outstream << "\n";
		missing = false;
		outstream << "INFORMATION ABOUT RESIDUE LEVEL TASK OPERATION \"" << component_name << "\":\n\n";
		utility::tag::XMLSchemaDefinition xsd;
		rltaskfactory->provide_xml_schema( component_name, xsd );
		outstream << xsd.human_readable_summary( component_name, "rlto" );
	}

	return missing;
}

/// @brief Prints out XSD information about the XML-accessible options for a given set of RosettaScipts-accessible
/// movers, filters, task operations, or residue selectors.
/// @details Calls the single string version.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
print_information(
	utility::vector1 < std::string > const &component_names
) {
	core::Size const ncomponents( component_names.size() );
	runtime_assert_string_msg( ncomponents > 0, "Error!  The rosetta_scripts application was used with the -parser:info flag, but nothing was provided after this flag.  The user must specify the name of at least one mover, filter, task operation, or residue selector for which to retrieve information." );

	utility::vector1 < std::string > failed_components;
	std::stringstream outstream("");
	outstream << "\nThe rosetta_scripts application was used with the -parser:info flag.\nWriting options for the indicated movers/filters/task operations/residue selectors:\n";

	for ( core::Size i(1); i<=ncomponents; ++i ) {
		outstream << "--------------------------------------------------------------------------------\n";
		if ( print_information( component_names[i], outstream ) ) {
			failed_components.push_back( component_names[i] );
			outstream << "\"" << component_names[i] << "\" not found!\n";
		}
	}
	outstream << "--------------------------------------------------------------------------------\n";
	if ( failed_components.size() > 0 ) {
		outstream << "Warning: the following are not movers, filters, task operations, or residue selectors; no information could be found for these:\n";
		for ( core::Size i(1), imax(failed_components.size()); i<=imax; ++i ) {
			outstream << failed_components[i] << "\n";
		}
	}
	outstream << "\nThe rosetta_scripts application will now exit.";
	TR << outstream.str() << std::endl;
	TR.flush();
}

/// @brief Saves the XSD to the given file.
void save_schema(  std::string const & filename ) {
	std::string schema;
	try {
		protocols::rosetta_scripts::RosettaScriptsParser parser;
		schema = parser.xsd_for_rosetta_scripts();
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		TR.Error << "An error was encountered while attempting to generate the XSD schema - it will not be saved.\n";
		TR.Error << e.msg() << "\n";
		return;
	}

	std::ofstream ofs( filename );
	ofs << schema << std::endl;
}

/// @details dock_design_scripting provides an xml-based scripting capability
/// to run rosetta movers and filters defined in a text file provided by the
/// user. A full documentation of dock_design_scripting is available at:
/// manual_doxygen/applications/app_dock_design.dox
int
main( int argc, char * argv [] )
{
	try{
		protocols::abinitio::ClassicAbinitio::register_options();
		// setup random numbers and options
		devel::init(argc, argv);
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( option[ parser::info ].user() ) { // If the -parser::info option is used, just print information about the requested component(s) and exit.
			print_information( option[ parser::info ]() );
		} else if ( option[ parser::output_schema ].user() ) {
			save_schema( option[ parser::output_schema ] );
		} else if ( ! option[ parser::protocol ].user() ) { // Just print a template script and exit if no input script is provided.
				print_template_script();
		} else { // If an input script has been provided, then we're not printing a template script and exiting.
			bool const view( option[ parser::view ] );
			protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.
			//That means that this can only be used with parser=true
			option[ jd2::dd_parser ].value( true ); // So here we fix that. jd2_parser app makes no sense without this option=true
			if ( !option[ jd2::ntrials ].user() ) {
				// when using rosetta_scripts we want ntrials to be set to 1 if the user forgot to specify. o/w infinite loops might
				// occur. We don't want ntrials to be set as default to 1, b/c other protocols might want it to behave differently
				option[ jd2::ntrials ].value( 1 );
			}

			if ( view ) {
				protocols::viewer::viewer_main( my_main );
			} else {
#ifdef BOINC
				protocols::jd2::BOINCJobDistributor::get_instance()->go( mover );
#else
#ifdef USEMPI
				protocols::jd2::MPIFileBufJobDistributor::get_instance()->go( mover );
#else
				protocols::jd2::JobDistributor::get_instance()->go( mover );
#endif
#endif
			}
		}
	} catch( utility::excn::EXCN_Base& excn ) {
		basic::Error()
			<< "ERROR: Exception caught by rosetta_scripts application:"
			<< excn << std::endl;
		std::exit( 1 );
	}
}

