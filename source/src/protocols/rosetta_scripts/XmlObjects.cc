// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rosetta_scripts/XmlObjects.cc
/// @brief Class to load objects from xml
/// @details There are three ways to use this class.
/// 1. You may pass an .xml file to the static method create_from_file which will return
///    an XmlObjects file where you can use the get* methods
/// 2. You may pass an xml script in text form to create_from_string which will return
///    an XmlObjects file where you can use the get* methods
/// 3. You may pass a single xml tag to the static methods static_get* which will return
///    an instance of whatever you instantiated
/// This class cannot be used to modify the behavior of the xml script contained in the mover
///   "ParserProtocol" (at present). If you wish to change that behavior, look into the
///   the script vars set with -parser::script_vars and used with %%script_var%%
/// @author Brian Coventry


// Unit headers
#include <protocols/rosetta_scripts/XmlObjects.hh>

// Package headers
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

// Project headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Numeric headers

//// C++ headers
#include <string>



// Construct tracer.
static basic::Tracer TR("protocols.rosetta_scripts.XmlObjects");

namespace protocols {
namespace rosetta_scripts {

using namespace core;


XmlObjects::XmlObjects()= default;


XmlObjects::~XmlObjects()= default;

XmlObjects::XmlObjects( XmlObjects const & src ) {
	*this = src;
}

XmlObjects&
XmlObjects::operator=( XmlObjects const & src ) {

	if ( this != & src ) {

		score_functions_ = src.score_functions_;
		residue_selectors_ = src.residue_selectors_;
		filters_ = src.filters_;
		movers_ = src.movers_;
		task_operations_ = src.task_operations_;

	}

	return *this;
}

/// @brief Parses an xml-formatted string and returns an XmlObjects
/// @details This function will wrap your script with <ROSETTASCRIPTS>
///   if needed and will add a <PROTOCOLS> section if needed. The
///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
XmlObjectsCOP
XmlObjects::create_from_string(
	std::string const & xml_text ) {
	pose::Pose pose;
	return create_from_string( xml_text, pose, basic::options::option );
}

/// @brief Parses an xml-formatted string and returns an XmlObjects
/// @details This function will wrap your script with <ROSETTASCRIPTS>
///   if needed and will add a <PROTOCOLS> section if needed. The
///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
XmlObjectsCOP
XmlObjects::create_from_string(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return create_from_string( xml_text, pose, basic::options::option );
}

/// @brief Parses an xml-formatted string and returns an XmlObjects
/// @details This function will wrap your script with <ROSETTASCRIPTS>
///   if needed and will add a <PROTOCOLS> section if needed. The
///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
XmlObjectsCOP
XmlObjects::create_from_string(
	std::string const & xml_text,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string prepared_text = xml_text;
	prepare_xml_text( prepared_text );

	RosettaScriptsParser parser;
	bool modified_pose = false;

	XmlObjectsOP objs( new XmlObjects() );

	ParsedProtocolOP parsed_protocol;
	// This try catch block should be removed once PyRosetta can correctly output the error
	try {
		parsed_protocol = parser.generate_mover_and_apply_to_pose_xml_string(
			pose, options, modified_pose, prepared_text, "", "", objs );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl;
		TR << "Error creating parsed protocol." << std::endl;
		throw e;
	}

	objs->movers_[ "ParsedProtocol" ] = parsed_protocol ;

	return objs;

}

/// @brief Parses an xml file and returns an XmlObjects
/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
XmlObjectsCOP
XmlObjects::create_from_file(
	std::string const & filename ) {
	pose::Pose pose;
	return create_from_file( filename, pose, basic::options::option );
}

/// @brief Parses an xml file and returns an XmlObjects
/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
XmlObjectsCOP
XmlObjects::create_from_file(
	std::string const & filename,
	core::pose::Pose & pose ) {
	return create_from_file( filename, pose, basic::options::option );
}

/// @brief Parses an xml file and returns an XmlObjects
/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
XmlObjectsCOP
XmlObjects::create_from_file(
	std::string const & filename,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	RosettaScriptsParser parser;

	bool modified_pose = false;

	XmlObjectsOP objs( new XmlObjects() );
	ParsedProtocolOP parsed_protocol;
	// This try catch block should be removed once PyRosetta can correctly output the error
	try {
		parsed_protocol = parser.generate_mover_and_apply_to_pose(
			pose, options, modified_pose, filename, "", "", objs );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl;
		TR << "Error creating parsed protocol." << std::endl;
		throw e;
	}

	objs->movers_[ "ParsedProtocol" ] = parsed_protocol;

	return objs;

}

/// @brief Wraps an xml script with <ROSETTASCRIPTS> and adds a <PROTOCOLS>
///   section if these don't exist already.
/// This function will fail if someone has a mover named ROSETTASCRIPTS or PROTOCOLS
///   or has those words in a comment. These string.find methods aren't the best
///   solution here, but without regular expressions (to catch things like < ROSETTASCRIPTS >),
///   correctly determing whether or not to intervene is complicated. -bcov
void
XmlObjects::prepare_xml_text( std::string & xml_text ) {
	if ( xml_text.find("ROSETTASCRIPTS") == std::string::npos ) {
		if ( xml_text.find("PROTOCOLS") == std::string::npos ) {
			xml_text += "\n<PROTOCOLS>\n</PROTOCOLS>\n";
		}
		xml_text = "<ROSETTASCRIPTS>\n" + xml_text + "</ROSETTASCRIPTS>\n";
	}
}


/// @brief Initialization function from the DataMaps
void
XmlObjects::init_from_maps(
	basic::datacache::DataMap & data,
	filters::Filters_map const & filters,
	moves::Movers_map const & movers ) {

	if ( data.has_type( "ResidueSelector" ) ) {
		for ( auto pair : data[ "ResidueSelector" ] ) {
			std::string name = pair.first;
			residue_selectors_[ name ] =
				data.get_ptr< select::residue_selector::ResidueSelector const >( "ResidueSelector", name );
		}
	}
	if ( data.has_type( "scorefxns" ) ) {
		for ( auto pair : data[ "scorefxns" ] ) {
			std::string name = pair.first;
			score_functions_[ name ] =
				data.get_ptr< scoring::ScoreFunction const >( "scorefxns", name );
		}
	}
	if ( data.has_type( "task_operations" ) ) {
		for ( auto pair : data[ "task_operations" ] ) {
			std::string name = pair.first;
			task_operations_[ name ] =
				data.get_ptr< pack::task::operation::TaskOperation >( "task_operations", name );
		}
	}

	for ( auto pair : filters ) {
		filters_[ pair.first ] = pair.second;
	}
	for ( auto pair : movers ) {
		movers_[ pair.first ] = pair.second;
	}

}

/// @brief Constructs a single ScoreFunction from xml
/// @details Pass this function a single <ScoreFunction /> tag and it will
///   return to you that ScoreFunction
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
core::scoring::ScoreFunctionOP
XmlObjects::static_get_score_function(
	std::string const & xml_text ) {
	pose::Pose pose;
	return static_get_score_function( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single ScoreFunction from xml
/// @details Pass this function a single <ScoreFunction /> tag and it will
///   return to you that ScoreFunction
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::scoring::ScoreFunctionOP
XmlObjects::static_get_score_function(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return static_get_score_function( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single ScoreFunction from xml
/// @details Pass this function a single <ScoreFunction /> tag and it will
///   return to you that ScoreFunction
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::scoring::ScoreFunctionOP
XmlObjects::static_get_score_function(
	std::string const & xml_text_in,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string xml_text = xml_text_in;
	std::string name;
	try {
		name = get_or_set_tag_name( xml_text );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl; // This line should be removed once PyRosetta correctly outputs the exception
		TR << "Error parsing tag. Remember to specify only one ScoreFunction!" << std::endl;
		throw e;
	}

	std::string wrapped_text = "<SCOREFXNS>" + xml_text + "</SCOREFXNS>";
	XmlObjectsCOP objs = create_from_string( wrapped_text, pose, options );

	return objs->get_score_function( name );
}

/// @brief Constructs a single ResidueSelector from xml
/// @details Pass this function a single ResidueSelector tag and it will
///   return to you that ResidueSelector. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
core::select::residue_selector::ResidueSelectorOP
XmlObjects::static_get_residue_selector(
	std::string const & xml_text ) {
	pose::Pose pose;
	return static_get_residue_selector( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single ResidueSelector from xml
/// @details Pass this function a single ResidueSelector tag and it will
///   return to you that ResidueSelector. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::select::residue_selector::ResidueSelectorOP
XmlObjects::static_get_residue_selector(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return static_get_residue_selector( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single ResidueSelector from xml
/// @details Pass this function a single ResidueSelector tag and it will
///   return to you that ResidueSelector. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::select::residue_selector::ResidueSelectorOP
XmlObjects::static_get_residue_selector(
	std::string const & xml_text_in,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string xml_text = xml_text_in;
	std::string name;
	try {
		name = get_or_set_tag_name( xml_text );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl; // This line should be removed once PyRosetta correctly outputs the exception
		TR << "Error parsing tag. Remember to specify only one ResidueSelector!" << std::endl;
		throw e;
	}

	std::string wrapped_text = "<RESIDUE_SELECTORS>" + xml_text + "</RESIDUE_SELECTORS>";
	XmlObjectsCOP objs = create_from_string( wrapped_text, pose, options );

	return objs->get_residue_selector( name );
}

/// @brief Constructs a single Filter from xml
/// @details Pass this function a single Filter tag and it will
///   return to you that Filter. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
protocols::filters::FilterOP
XmlObjects::static_get_filter(
	std::string const & xml_text ) {
	pose::Pose pose;
	return static_get_filter( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single Filter from xml
/// @details Pass this function a single Filter tag and it will
///   return to you that Filter. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
protocols::filters::FilterOP
XmlObjects::static_get_filter(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return static_get_filter( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single Filter from xml
/// @details Pass this function a single Filter tag and it will
///   return to you that Filter. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
protocols::filters::FilterOP
XmlObjects::static_get_filter(
	std::string const & xml_text_in,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string xml_text = xml_text_in;
	std::string name;
	try {
		name = get_or_set_tag_name( xml_text );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl; // This line should be removed once PyRosetta correctly outputs the exception
		TR << "Error parsing tag. Remember to specify only one Filter!" << std::endl;
		throw e;
	}

	std::string wrapped_text = "<FILTERS>" + xml_text + "</FILTERS>";
	XmlObjectsCOP objs = create_from_string( wrapped_text, pose, options );

	return objs->get_filter( name );
}

/// @brief Constructs a single Mover from xml
/// @details Pass this function a single Mover tag and it will
///   return to you that Mover. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
protocols::moves::MoverOP
XmlObjects::static_get_mover(
	std::string const & xml_text ) {
	pose::Pose pose;
	return static_get_mover( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single Mover from xml
/// @details Pass this function a single Mover tag and it will
///   return to you that Mover. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
protocols::moves::MoverOP
XmlObjects::static_get_mover(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return static_get_mover( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single Mover from xml
/// @details Pass this function a single Mover tag and it will
///   return to you that Mover. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
protocols::moves::MoverOP
XmlObjects::static_get_mover(
	std::string const & xml_text_in,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string xml_text = xml_text_in;
	std::string name;
	try {
		name = get_or_set_tag_name( xml_text );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl; // This line should be removed once PyRosetta correctly outputs the exception
		TR << "Error parsing tag. Remember to specify only one Mover!" << std::endl;
		throw e;
	}

	std::string wrapped_text = "<MOVERS>" + xml_text + "</MOVERS>";
	XmlObjectsCOP objs = create_from_string( wrapped_text, pose, options );

	return objs->get_mover( name );
}

/// @brief Constructs a single TaskOperation from xml
/// @details Pass this function a single TaskOperation tag and it will
///   return to you that TaskOperation. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   A pose is needed for the initialization of some Movers and may
///   be required if this functions throws an error.
core::pack::task::operation::TaskOperationOP
XmlObjects::static_get_task_operation(
	std::string const & xml_text ) {
	pose::Pose pose;
	return static_get_task_operation( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single TaskOperation from xml
/// @details Pass this function a single TaskOperation tag and it will
///   return to you that TaskOperation. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::pack::task::operation::TaskOperationOP
XmlObjects::static_get_task_operation(
	std::string const & xml_text,
	core::pose::Pose & pose ) {
	return static_get_task_operation( xml_text, pose, basic::options::option );
}

/// @brief Constructs a single TaskOperation from xml
/// @details Pass this function a single TaskOperation tag and it will
///   return to you that TaskOperation. For C++ users, you may need
///   to use std::dynamic_pointer_cast<   > after this call.
///   The pose is needed for the initialization of some Movers and may
///   be modified if an APPLY_TO_POSE section is present.
core::pack::task::operation::TaskOperationOP
XmlObjects::static_get_task_operation(
	std::string const & xml_text_in,
	core::pose::Pose & pose,
	utility::options::OptionCollection const & options ) {

	std::string xml_text = xml_text_in;
	std::string name;
	try {
		name = get_or_set_tag_name( xml_text );
	} catch ( utility::excn::Exception const & e ) {
		TR << e.msg() << std::endl; // This line should be removed once PyRosetta correctly outputs the exception
		TR << "Error parsing tag. Remember to specify only one TaskOperation!" << std::endl;
		throw e;
	}

	std::string wrapped_text = "<TASKOPERATIONS>" + xml_text + "</TASKOPERATIONS>";
	XmlObjectsCOP objs = create_from_string( wrapped_text, pose, options );

	return objs->get_task_operation( name );
}

/// @brief List all the ScoreFunctions contained by name
utility::vector1< std::string >
XmlObjects::list_score_functions() const {
	utility::vector1< std::string > names;
	for ( auto pair : score_functions_ ) {
		names.push_back( pair.first );
	}
	return names;
}

/// @brief List all the ResidueSelectors contained by name
utility::vector1< std::string >
XmlObjects::list_residue_selectors() const {
	utility::vector1< std::string > names;
	for ( auto pair : residue_selectors_ ) {
		names.push_back( pair.first );
	}
	return names;
}

/// @brief List all the Filters contained by name
utility::vector1< std::string >
XmlObjects::list_filters() const {
	utility::vector1< std::string > names;
	for ( auto pair : filters_ ) {
		names.push_back( pair.first );
	}
	return names;
}

/// @brief List all the Movers contained by name
utility::vector1< std::string >
XmlObjects::list_movers() const {
	utility::vector1< std::string > names;
	for ( auto pair : movers_ ) {
		names.push_back( pair.first );
	}
	return names;
}

/// @brief List all the TaskOperations by name
utility::vector1< std::string >
XmlObjects::list_task_operations() const {
	utility::vector1< std::string > names;
	for ( auto pair : task_operations_ ) {
		names.push_back( pair.first );
	}
	return names;
}


/// @brief Extract a ScoreFunction by name after using one of the create* methods
/// @details This returns a clone and cannot be used to modify the ParsedProtocol
scoring::ScoreFunctionOP
XmlObjects::get_score_function( std::string const & name ) const {
	if ( score_functions_.count( name ) == 0 ) {
		throw std::runtime_error( "ScoreFunction with name \"" + name +
			"\" does not exist!!" );
	}
	// DO NOT switch this function to return a COP and remove this clone. PyRosetta
	// allows one to modify COPs which allows one to change the const private data of this
	// class and could even allow one to modify the ParsedProtocol. This behavior is
	// side-effect programming and must be disallowed
	return score_functions_.at( name )->clone();
}

/// @brief Extract a ResidueSelector by name after using one of the create* methods
/// @details This returns a clone and cannot be used to modify the ParsedProtocol
select::residue_selector::ResidueSelectorOP
XmlObjects::get_residue_selector( std::string const & name ) const {
	if ( residue_selectors_.count( name ) == 0 ) {
		throw std::runtime_error( "ResidueSelector with name \"" + name +
			"\" does not exist!!" );
	}
	// DO NOT switch this function to return a COP and remove this clone. PyRosetta
	// allows one to modify COPs which allows one to change the const private data of this
	// class and could even allow one to modify the ParsedProtocol. This behavior is
	// side-effect programming and must be disallowed
	return residue_selectors_.at( name )->clone();
}


/// @brief Extract a Filter by name after using one of the create* methods
/// @details This returns a clone and cannot be used to modify the ParsedProtocol
filters::FilterOP
XmlObjects::get_filter( std::string const & name ) const {
	if ( filters_.count( name ) == 0 ) {
		throw std::runtime_error( "Filter with name \"" + name +
			"\" does not exist!!" );
	}
	// DO NOT switch this function to return a COP and remove this clone. PyRosetta
	// allows one to modify COPs which allows one to change the const private data of this
	// class and could even allow one to modify the ParsedProtocol. This behavior is
	// side-effect programming and must be disallowed
	return filters_.at( name )->clone();
}


/// @brief Extract a Mover by name after using one of the create* methods
/// @details This returns a clone and cannot be used to modify the ParsedProtocol
moves::MoverOP
XmlObjects::get_mover( std::string const & name ) const {
	if ( movers_.count( name ) == 0 ) {
		throw std::runtime_error( "Mover with name \"" + name +
			"\" does not exist!!" );
	}
	// DO NOT switch this function to return a COP and remove this clone. PyRosetta
	// allows one to modify COPs which allows one to change the const private data of this
	// class and could even allow one to modify the ParsedProtocol. This behavior is
	// side-effect programming and must be disallowed
	return movers_.at( name )->clone();
}


/// @brief Extract a TaskOperation by name after using one of the create* methods
/// @details This returns a clone and cannot be used to modify the ParsedProtocol
pack::task::operation::TaskOperationOP
XmlObjects::get_task_operation( std::string const & name ) const {
	if ( task_operations_.count( name ) == 0 ) {
		throw std::runtime_error( "TaskOperation with name \"" + name +
			"\" does not exist!!" );
	}
	// DO NOT switch this function to return a COP and remove this clone. PyRosetta
	// allows one to modify COPs which allows one to change the const private data of this
	// class and could even allow one to modify the ParsedProtocol. This behavior is
	// side-effect programming and must be disallowed
	return task_operations_.at( name )->clone();
}

std::string
XmlObjects::get_or_set_tag_name( std::string & tag_string ) {
	utility::tag::TagOP tag = utility::tag::Tag::create( tag_string );
	std::string name;
	if ( tag->hasOption( "name" ) ) {
		name = tag->getOption<std::string>( "name", "" );
	} else {
		name = "XmlObject_name";
		tag->setOption( "name", name );
		std::stringstream ss;
		tag->write( ss );
		tag_string = ss.str();
	}
	return name;
}


/// @brief Generate string representation of XmlObjects for debugging purposes.
void
XmlObjects::show( std::ostream & output ) const {

	output << "ScoreFunctions:" << std::endl;
	for ( auto pair : score_functions_ ) {
		output << "  " << pair.first << std::endl;
	}
	output << "ResidueSelectors:" << std::endl;
	for ( auto pair : residue_selectors_ ) {
		output << "  " << pair.first << std::endl;
	}
	output << "TaskOperations:" << std::endl;
	for ( auto pair : task_operations_ ) {
		output << "  " << pair.first << std::endl;
	}
	output << "Filters:" << std::endl;
	for ( auto pair : filters_ ) {
		output << "  " << pair.first << std::endl;
	}
	output << "Movers:" << std::endl;
	for ( auto pair : movers_ ) {
		output << "  " << pair.first << std::endl;
	}
}

/// @brief Insertion operator (overloaded so that XmlObjects can be "printed" in PyRosetta).
std::ostream &
operator << ( std::ostream & output, XmlObjects const & object_to_output ) {
	object_to_output.show( output );
	return output;
}



} // rosetta_scripts
} // protocols


