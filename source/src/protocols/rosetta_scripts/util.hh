// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/RosettaScripts/util.hh
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
///         Jacob Corn (jecorn@u.washington.edu)
///         Rocco Moretti (rmoretti@u.washington.edu)
///         Eva-Maria Strauch (evas01@uw.edu)


#ifndef INCLUDED_protocols_rosetta_scripts_util_hh
#define INCLUDED_protocols_rosetta_scripts_util_hh

// Unit headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>

// Utillity Headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace rosetta_scripts {


/// @brief find source residue that is nearest to res on source. If distance is greater than 2.0A, return 0. chain=0, search all chains, chain=1,2,3 etc. search only that chain
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain = 0 );

/// @brief find nearest residue and also return the distance
void
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size & target_res, core::Real & dist, core::Size const chain = 0 );

/// @brief finds the nearest disulife to given residue on pose by searching both up and down stream to the given postion
core::Size
find_nearest_disulfide( core::pose::Pose const & pose, core::Size const res);

/// @brief returns a vector containing all the residues with a given packer state according to the TF
utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/ );



///////////////////////////////////////////////////////////
//////////////////// Task Operations //////////////////////

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

/// allows the transfer of whole taskfactories on the datamap. This way a "base" taskfactory can be created, transferred on the datamap, and
/// individual mover's specific taskoperations can be added on top
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap /*const*/ & data, core::pack::task::TaskFactoryOP & task_factory );

core::pack::task::TaskFactoryOP
parse_task_operations( std::string const & task_list, basic::datacache::DataMap const & data );


///////////////////// Attributes /////////////////////////

/// @brief Appends the 'task_operation' attribute
/// @details "description" can be used to specify for what the TaskOperations are being used for.
void
attributes_for_parse_task_operations(
	utility::tag::AttributeList & attributes,
	std::string const & description = "" );

/// @brief Appends the 'task_operation' and 'task_factory' attributes.
void
attributes_for_parse_task_operations_w_factory(
	utility::tag::AttributeList & attributes,
	std::string const & used_for_descr = "" );


/////////////////////////////////////////////////////////////
//////////////////// Residue Selectors //////////////////////

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" option in tag
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
/// @note The default option is "residue_selector".  However, this function can be used to
/// get selectors with other option names if another string is passed for the third parameter.
/// This is useful in cases with multiple selectors (e.g. "first_selector", "second_selector",
/// etc.).
core::select::residue_selector::ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name = "residue_selector" );

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data );


///////////////////// Attributes ///////////////////////////

/// @brief Appends the attributes read by parse_residue_selector
void
attributes_for_parse_residue_selector( utility::tag::AttributeList & attributes, std::string const & description = "" );


///////////////////////////////////////////////////////////
//////////////////// ScoreFunction ////////////////////////

/// @brief Look up the score function defined in the <SCOREFXNS/>
/// through the given option. Defaults to 'commandline'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="commandline" );

/// @brief Look up the score function defined in the <SCOREFXNS/>
///through the option 'scorefxn='. Defaults to 'commandline'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="commandline" );

/// @brief Look up the name of assigned score function to the given
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name);

/// @brief Look up the name of assigned score function to the 'scorefxn='
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag);


///////////////////// Attributes ///////////////////////////

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes );

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes,
	std::string const & option_name);

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description );

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & option_name,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes);

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name );

/// @brief Appends the attributes read by parse_score_function with description
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function w/ name argument and description
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function w/ name argument and description.
/// @details This version appends the attributes as required attributes.
/// @author Vikram K. Mulligan.
void
attributes_for_parse_score_function_w_description_when_required( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description = ""
);


/////////////////////////////////////////////////////////
//////////////////// MoveMap ////////////////////////////

/// @brief Parse a MoveMap factory from tags using the old MoveMap specification syntax.
/// @details Will return nullptr if the tag doesn't contain a MoveMap specification.
core::select::movemap::MoveMapFactoryOP
parse_movemap_factory_legacy(
	utility::tag::TagCOP in_tag,
	basic::datacache::DataMap & data,
	bool const reset_movemap = true, /* should we turn everything to true at start?*/
	core::select::movemap::MoveMapFactoryOP mmf_to_modify = nullptr );


/// @brief Parse a MoveMap factory from tags using the old MoveMap specification syntax.
/// Will look for a MoveMap entry with a particular name entry.
/// @details Will return nullptr if the tag doesn't contain a MoveMap specification.
core::select::movemap::MoveMapFactoryOP
parse_named_movemap_factory_legacy(
	utility::tag::TagCOP in_tag,
	std::string const & name,
	basic::datacache::DataMap & data,
	bool const reset_movemap = true, // should we turn everything to true at start?
	core::select::movemap::MoveMapFactoryOP mmf_to_modify = nullptr );

///////////////////// Attributes ///////////////////////////

/// @brief Adds a subelement to an input subelement list describing a MoveMap subtag
/// that will be used by the parse_movemap_legacy function.
void
append_subelement_for_parse_movemap_factory_legacy(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::XMLSchemaSimpleSubelementList & subelements,
	std::string const & description = ""
);

///////////////////////////////////////////////////////////
//////////////////// Reference Pose ///////////////////////

core::pose::PoseOP
saved_reference_pose(
	utility::tag::TagCOP in_tag,
	basic::datacache::DataMap & data_map,
	std::string const & tag_str="reference_name" );

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
/// Will look for tag in in_tag variable
void
attributes_for_saved_reference_pose( utility::tag::AttributeList & attributes, std::string const & attribute_name="reference_name" );

/// @brief convenience function to access pointers to poses that will be stored
/// in the data map at an arbitrary point during an RS protocol
/// Will look for tag in in_tag variable
void
attributes_for_saved_reference_pose_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description,
	std::string const & attribute_name="reference_name" );

/////////////////////////////////////////////////////////
//////////////////// Filter ////////////////////////////

protocols::filters::FilterOP
parse_filter( std::string const & filter_name, protocols::filters::Filters_map const & d );


/////////////////////////////////////////////////////////
//////////////////// Mover //////////////////////////////
protocols::moves::MoverOP
parse_mover( std::string const & mover_name, protocols::moves::Movers_map const & d );


/////////////////////////////////////////////////////////
//////////////////// XYZ Vector /////////////////////////
numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagCOP xyz_vector_tag );

void
attributes_for_parse_xyz_vector( utility::tag::AttributeList & attlist );

/////////////////////////////////////////////////////////
//////////////DATABASE SESSIONS//////////////////////////

utility::sql_database::sessionOP
parse_database_session(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & datamap
);

void
attributes_for_parse_database_session(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::AttributeList & attlist
);


///This is kind of a strange place for this, but for library-level reasons it needs to be more accessible than a more logical home with ReportToDB, and cannot live in basic because it needs other functions in this file.  (There is also value in not creating a new file b/c it breaks the fast-compile system XML XSD XRW is using, and it's 6pm on Friday!)

//void
//attributes_for_report_to_db( utility::tag::AttributeList &, utility::tag::XMLSchemaDefinition & );




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



} // RosettaScripts
} // protocols

#endif /*INCLUDED_protocols_RosettaScripts_util_HH*/
