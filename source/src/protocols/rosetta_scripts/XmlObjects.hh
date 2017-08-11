// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rosetta_scripts/XmlObjects.hh
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


#ifndef INCLUDED_protocols_rosetta_scripts_XmlObjects_hh
#define INCLUDED_protocols_rosetta_scripts_XmlObjects_hh

// Unit headers
#include <protocols/rosetta_scripts/XmlObjects.fwd.hh>

// Project headers
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Numeric headers


// C++ headers
#include <iostream>
#include <string>


namespace protocols {
namespace rosetta_scripts {


class XmlObjects : public utility::pointer::ReferenceCount {
public:

	XmlObjects();
	~XmlObjects() override;

	XmlObjects( XmlObjects const & src );
	XmlObjects& operator=( XmlObjects const & src );

	/// @brief Parses an xml-formatted string and returns an XmlObjects
	/// @details This function will wrap your script with <ROSETTASCRIPTS>
	///   if needed and will add a <PROTOCOLS> section if needed. The
	///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	XmlObjectsCOP
	create_from_string(
		std::string const & xml_text );

	/// @brief Parses an xml-formatted string and returns an XmlObjects
	/// @details This function will wrap your script with <ROSETTASCRIPTS>
	///   if needed and will add a <PROTOCOLS> section if needed. The
	///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	XmlObjectsCOP
	create_from_string(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Parses an xml-formatted string and returns an XmlObjects
	/// @details This function will wrap your script with <ROSETTASCRIPTS>
	///   if needed and will add a <PROTOCOLS> section if needed. The
	///   ParsedProtocol (xml script mover) is available as "ParsedProtocol".
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	XmlObjectsCOP
	create_from_string(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Parses an xml file and returns an XmlObjects
	/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	XmlObjectsCOP
	create_from_file(
		std::string const & filename );

	/// @brief Parses an xml file and returns an XmlObjects
	/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	XmlObjectsCOP
	create_from_file(
		std::string const & filename,
		core::pose::Pose & pose );

	/// @brief Parses an xml file and returns an XmlObjects
	/// @details The ParsedProtocol (xml script mover) is available as "ParsedProtocol"
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	XmlObjectsCOP
	create_from_file(
		std::string const & filename,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Initialization function from the DataMaps
	void
	init_from_maps(
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers );

	/// @brief Constructs a single ScoreFunction from xml
	/// @details Pass this function a single <ScoreFunction /> tag and it will
	///   return to you that ScoreFunction
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::scoring::ScoreFunctionOP
	static_get_score_function(
		std::string const & xml_text );

	/// @brief Constructs a single ScoreFunction from xml
	/// @details Pass this function a single <ScoreFunction /> tag and it will
	///   return to you that ScoreFunction
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::scoring::ScoreFunctionOP
	static_get_score_function(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Constructs a single ScoreFunction from xml
	/// @details Pass this function a single <ScoreFunction /> tag and it will
	///   return to you that ScoreFunction
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::scoring::ScoreFunctionOP
	static_get_score_function(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Constructs a single ResidueSelector from xml
	/// @details Pass this function a single ResidueSelector tag and it will
	///   return to you that ResidueSelector. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	core::select::residue_selector::ResidueSelectorOP
	static_get_residue_selector(
		std::string const & xml_text );

	/// @brief Constructs a single ResidueSelector from xml
	/// @details Pass this function a single ResidueSelector tag and it will
	///   return to you that ResidueSelector. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::select::residue_selector::ResidueSelectorOP
	static_get_residue_selector(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Constructs a single ResidueSelector from xml
	/// @details Pass this function a single ResidueSelector tag and it will
	///   return to you that ResidueSelector. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::select::residue_selector::ResidueSelectorOP
	static_get_residue_selector(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Constructs a single Filter from xml
	/// @details Pass this function a single Filter tag and it will
	///   return to you that Filter. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	protocols::filters::FilterOP
	static_get_filter(
		std::string const & xml_text );

	/// @brief Constructs a single Filter from xml
	/// @details Pass this function a single Filter tag and it will
	///   return to you that Filter. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	protocols::filters::FilterOP
	static_get_filter(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Constructs a single Filter from xml
	/// @details Pass this function a single Filter tag and it will
	///   return to you that Filter. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	protocols::filters::FilterOP
	static_get_filter(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Constructs a single Mover from xml
	/// @details Pass this function a single Mover tag and it will
	///   return to you that Mover. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	protocols::moves::MoverOP
	static_get_mover(
		std::string const & xml_text );

	/// @brief Constructs a single Mover from xml
	/// @details Pass this function a single Mover tag and it will
	///   return to you that Mover. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	protocols::moves::MoverOP
	static_get_mover(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Constructs a single Mover from xml
	/// @details Pass this function a single Mover tag and it will
	///   return to you that Mover. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	protocols::moves::MoverOP
	static_get_mover(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief Constructs a single TaskOperation from xml
	/// @details Pass this function a single TaskOperation tag and it will
	///   return to you that TaskOperation. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   A pose is needed for the initialization of some Movers and may
	///   be required if this functions throws an error.
	static
	core::pack::task::operation::TaskOperationOP
	static_get_task_operation(
		std::string const & xml_text );

	/// @brief Constructs a single TaskOperation from xml
	/// @details Pass this function a single TaskOperation tag and it will
	///   return to you that TaskOperation. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::pack::task::operation::TaskOperationOP
	static_get_task_operation(
		std::string const & xml_text,
		core::pose::Pose & pose );

	/// @brief Constructs a single TaskOperation from xml
	/// @details Pass this function a single TaskOperation tag and it will
	///   return to you that TaskOperation. For C++ users, you may need
	///   to use std::dynamic_pointer_cast<   > after this call.
	///   The pose is needed for the initialization of some Movers and may
	///   be modified if an APPLY_TO_POSE section is present.
	static
	core::pack::task::operation::TaskOperationOP
	static_get_task_operation(
		std::string const & xml_text,
		core::pose::Pose & pose,
		utility::options::OptionCollection const & options );

	/// @brief List all the ScoreFunctions contained by name
	utility::vector1< std::string >
	list_score_functions() const;

	/// @brief List all the ResidueSelectors contained by name
	utility::vector1< std::string >
	list_residue_selectors() const;

	/// @brief List all the Filters contained by name
	utility::vector1< std::string >
	list_filters() const;

	/// @brief List all the Movers contained by name
	utility::vector1< std::string >
	list_movers() const;

	/// @brief List all the TaskOperations by name
	utility::vector1< std::string >
	list_task_operations() const;

	/// @brief Extract a ScoreFunction by name after using one of the create* methods
	/// @details This returns a clone and cannot be used to modify the ParsedProtocol
	core::scoring::ScoreFunctionOP
	get_score_function( std::string const & name ) const;

	/// @brief Extract a ResidueSelector by name after using one of the create* methods
	/// @details This returns a clone and cannot be used to modify the ParsedProtocol
	core::select::residue_selector::ResidueSelectorOP
	get_residue_selector( std::string const & name ) const;

	/// @brief Extract a Filter by name after using one of the create* methods
	/// @details This returns a clone and cannot be used to modify the ParsedProtocol
	protocols::filters::FilterOP
	get_filter( std::string const & name ) const;

	/// @brief Extract a Mover by name after using one of the create* methods
	/// @details This returns a clone and cannot be used to modify the ParsedProtocol
	protocols::moves::MoverOP
	get_mover( std::string const & name ) const;

	/// @brief Extract a TaskOperation by name after using one of the create* methods
	/// @details This returns a clone and cannot be used to modify the ParsedProtocol
	core::pack::task::operation::TaskOperationOP
	get_task_operation( std::string const & name ) const;


	/// @brief  Generate string representation of XmlObjects for debugging purposes.
	void show( std::ostream & output = std::cout ) const;

	// Insertion operator (overloaded so that XmlObjects can be "printed" in PyRosetta).
	friend std::ostream & operator << ( std::ostream & output, XmlObjects const & object_to_output );

private:

	/// @brief Gets a tag's name field or sets it to XmlObject_name
	static
	std::string
	get_or_set_tag_name( std::string & tag_string );

	/// @brief Wraps an xml script with <ROSETTASCRIPTS> and adds a <PROTOCOLS>
	///   section if these don't exist already.
	static
	void
	prepare_xml_text( std::string & xml_text );


	std::map< std::string, core::scoring::ScoreFunctionCOP > score_functions_;
	std::map< std::string, core::select::residue_selector::ResidueSelectorCOP > residue_selectors_;
	std::map< std::string, protocols::filters::FilterCOP > filters_;
	std::map< std::string, protocols::moves::MoverCOP > movers_;
	std::map< std::string, core::pack::task::operation::TaskOperationCOP > task_operations_;


};

} // rosetta_scripts
} // protocols


#endif
