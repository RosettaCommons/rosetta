// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/RosettaScripts/util.cc
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu),
///     Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/task/PackerTask.hh>

// Project Headers
#include <core/types.hh>
#include <core/init/score_function_corrections.hh>
#include <protocols/filters/Filter.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <core/id/types.hh>
#include <core/chemical/VariantType.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>
#include <core/select/jump_selector/JumpIndexSelector.hh>
#include <core/select/residue_selector/ResidueSpanSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>

#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/sql_database/types.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

#include <fstream>

static basic::Tracer TR( "protocols.RosettaScripts.util" );

namespace protocols {
namespace rosetta_scripts {

using namespace core::scoring;
using namespace protocols::moves;
using namespace core;
using namespace std;
using utility::vector1;


/// @brief Return the number of the residue on source that is nearest to res on target. If the distance
/// is greater than 2.0 returns 0 to indicate error
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ ){
	//TR<<"looking for neiboughrs of: "<<source.pdb_info()->name()<< " and residue "<<res<<std::endl;
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	for ( core::Size i = 1; i <= source.size(); ++i ) {
		if ( source.residue( i ).is_ligand() ) continue;
		if ( chain && source.residue( i ).chain() != chain ) continue;
		core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
		if ( dist <= min_dist ) {
			min_dist = dist;
			nearest_res = i;
		}
	}
	if ( min_dist <= 3.0 ) return nearest_res;
	else return 0;
}

void
find_nearest_res(  core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size & target_res, core::Real & target_dist, core::Size const chain/*=0*/ ){
	target_res = 0; target_dist = 0.0;
	core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
	for ( core::Size i = 1; i <= source.size(); ++i ) {
		if ( source.residue( i ).is_ligand() ) continue;
		if ( chain && source.residue( i ).chain() != chain ) continue;
		core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
		if ( dist <= min_dist ) {
			min_dist = dist;
			nearest_res = i;
		}
	}
	if ( min_dist <= 3.0 ) {
		target_res = nearest_res;
		target_dist = min_dist;
	}
}


utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/) {
	utility::vector1< core::Size > designable_vec, packable_vec, both;
	designable_vec.clear(); packable_vec.clear(); both.clear();
	core::pack::task::PackerTaskOP packer_task( tf->create_task_and_apply_taskoperations( pose ) );
	for ( core::Size resi=1; resi<=pose.size(); ++resi ) {
		if ( packer_task->being_designed( resi ) ) {
			designable_vec.push_back( resi );
		} else if ( packer_task->being_packed( resi ) ) {
			packable_vec.push_back( resi );
		}
	}
	if ( designable && packable ) {
		both.insert( both.begin(), designable_vec.begin(), designable_vec.end() );
		both.insert( both.end(), packable_vec.begin(), packable_vec.end() );
		return both;
	}
	if ( designable ) {
		return designable_vec;
	}
	return packable_vec;
}
/// @brief finds the nearest disulife to given residue on pose by searching both up and down stream to the given postion
core::Size
find_nearest_disulfide( core::pose::Pose const & pose, core::Size const res)
{
	core::Size disulfideN=0,disulfideC=0;
	for ( core::Size i = res; i <= pose.size(); ++i ) {
		if ( pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
			disulfideC=i;
			break;
		}
	}
	// TR<<"C-ter disulfide: "<<disulfideC<<std::endl;
	for ( core::Size i = res ; i > 0; --i ) {
		if ( pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
			disulfideN=i;
			break;
		}
	}
	//TR<<"N-ter disulfide: "<<disulfideN<<std::endl;
	if ( (disulfideN==0)&&(disulfideC==0) ) {
		utility_exit_with_message("Could not find disulfides on: "+pose.pdb_info()->name());
	}
	if ( ((disulfideC-res)>(res-disulfideN))&&disulfideN!=0 ) {
		return disulfideN;
	}

	return disulfideC;
}


///////////////////////////////////////////////////////////
//////////////////// Task Operations //////////////////////


/// @details This is essentially a shameless copy of Justin's PackRotamersMover::parse_task_operations. In truth
/// DesignRepackMover should disappear into Justin's better organized class, but this will wait... (SJF)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP new_task_factory( new TaskFactory );
	if ( ! tag->hasOption("task_operations") ) return new_task_factory;
	TR<<"Object "<<tag->getOption< std::string >( "name", "no_name" )<<" reading the following task_operations: ";
	return( parse_task_operations( tag->getOption< std::string >( "task_operations"), data ) );
}

core::pack::task::TaskFactoryOP
parse_task_operations( std::string const & task_list, basic::datacache::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP new_task_factory( new TaskFactory );
	std::string const t_o_val( task_list );
	using StringVec = utility::vector1<std::string>;
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations\n";
	for ( auto const & t_o_key : t_o_keys ) {
		if ( data.has( "task_operations", t_o_key ) ) {
			new_task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR<<t_o_key<<' ';
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	TR<<std::endl;
	return new_task_factory;
}

///option to add or refer to a Taskfactory through the datamap, similar to how to add/refer to movemap OPs (EMS)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap & data, core::pack::task::TaskFactoryOP & task_factory /*, bool const reset_taskfactory */)
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	if ( tag->hasOption("task_factory" ) ) {
		std::string const name( tag->getOption< std::string >("task_factory") );
		TR <<"taskfacotry name: " << name << std::endl;

		if ( data.has( "TaskFactory", name ) ) {
			task_factory = data.get_ptr<TaskFactory>( "TaskFactory", name );
			TR<<"found helper task_factory: "<< name <<" for mover: "<<tag->getName()<< std::endl;
		} else { // ( !data.has( "TaskFactory", name ) ){
			std::string tf_string = "TaskFactory";
			task_factory = core::pack::task::TaskFactoryOP( new TaskFactory );
			data.add( tf_string , name , task_factory );
			TR<<"adding new TaskFactory to the datamap: "<< name  <<std::endl;
		}
	}

	if ( ! tag->hasOption("task_operations") ) return task_factory;

	std::string const t_o_val( tag->getOption<std::string>("task_operations") );
	using StringVec = utility::vector1<std::string>;
	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations to mover "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";

	for ( auto const & t_o_key : t_o_keys ) {

		if ( data.has( "task_operations", t_o_key ) ) {
			task_factory->push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			TR<<t_o_key<<' ';
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
		}
	}
	TR<<std::endl;
	return task_factory;
}

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	using core::pack::task::operation::TaskOperationOP;
	using core::pack::task::operation::TaskOperation;
	using String = std::string;

	utility::vector1< TaskOperationOP > task_operations;
	String const t_o_val( tag->getOption<String>("task_operations", "" ) );
	if ( t_o_val != "" ) {
		utility::vector1< String > const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( auto const & t_o_key : t_o_keys ) {
			if ( data.has( "task_operations", t_o_key ) ) {
				task_operations.push_back( data.get_ptr< TaskOperation >( "task_operations", t_o_key ) );
			} else {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
			}
		}
	}
	return task_operations;
}

////////////////////////////////////////////////////////////////////
///  Get attributes ( i.e. options ) for movers to build xml schemas

/// @brief Appends the attributes read by parse_task_operations
void
attributes_for_parse_task_operations( utility::tag::AttributeList & attributes, std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute( "task_operations", utility::tag::xsct_task_operation_comma_separated_list,
		( description == ""  ? "A comma separated list of TaskOperations to use." : description )  );
}


/// @brief Appends the attributes read by parse_task_operations when handed a TaskFactory
void
attributes_for_parse_task_operations_w_factory( utility::tag::AttributeList & attributes, std::string const & used_for_descr )
{
	std::string to_descrip = "A comma separated list of TaskOperations to use.";
	std::string tf_descrip = "A TaskFactory specification to use";
	if ( ! used_for_descr.empty() ) {
		to_descrip += " for " + used_for_descr;
		tf_descrip += " for " + used_for_descr;
	}

	attributes
		+ utility::tag::XMLSchemaAttribute( "task_operations", utility::tag::xsct_task_operation_comma_separated_list, to_descrip )
		+ utility::tag::XMLSchemaAttribute( "task_factory", utility::tag::xs_string , tf_descrip );
}

/////////////////////////////////////////////////////////////
//////////////////// Residue Selectors //////////////////////

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" option in tag
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
core::select::residue_selector::ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name )
{
	return core::select::residue_selector::parse_residue_selector( tag, data, option_name );
}

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data )
{
	return core::select::residue_selector::get_residue_selector( selector_name, data );
}


/// @brief Appends the attributes read by parse_residue_selector
void
attributes_for_parse_residue_selector( utility::tag::AttributeList & attributes, std::string const & description )
{
	return core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name(attributes, description );
}

///////////////////////////////////////////////////////////
//////////////////// ScoreFunction ////////////////////////

/// @details Utility function to find a scorefunction from parser-provided data.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	std::string const scorefxn_key( tag->getOption<std::string>(option_name, dflt_key) );
	if ( ! data.has( "scorefxns", scorefxn_key ) ) {
		if ( ! core::init::check_score_function_sanity( basic::options::option, scorefxn_key ) ) {
			// Incompatible default score functions aren't loaded.
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " was requested with incompatible options.");
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap. To add a score function to the data map, define a score function in the <SCOREFXNS/>.'");
		}
	}

	// We don't check "sanity" for the scorefxn_key, as that's user defined, and not necessarily representative of the actual scorefunction.

	return data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_key );
}



/// @details Utility function to find a scorefunction from
/// parser-provided data for the option 'scorefxn'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	return parse_score_function(tag, "scorefxn", data, dflt_key);
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name
) {
	return tag->getOption<std::string>(option_name, "commandline");
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag
) {
	return get_score_function_name(tag, "scorefxn");
}

///  Get attributes ( i.e. options ) for movers to build xml schemas
void
attributes_for_get_score_function_name( utility::tag::AttributeList & attributes )
{
	attributes_for_get_score_function_name_w_description( attributes, "" );
}

void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes,
	std::string const & option_name )
{
	attributes_for_get_score_function_name_w_description( attributes, option_name, "" );
}

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_get_score_function_name_w_description( attributes, "scorefxn", description);
}

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute( option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes )
{
	attributes_for_parse_score_function_w_description( attributes, "" );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes, std::string const & sfxn_option_name )
{
	attributes_for_parse_score_function_w_description( attributes, sfxn_option_name, "" );
}

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_parse_score_function_w_description( attributes, "scorefxn", description );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute(
		sfxn_option_name, utility::tag::xs_string ,
		( description == "" ? "Name of score function to use" : description ) );
}

/// @brief Appends the attributes read by parse_score_function w/ name argument and description.
/// @details This version appends the attributes as required attributes.
/// @author Vikram K. Mulligan.
void
attributes_for_parse_score_function_w_description_when_required( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description
) {
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		sfxn_option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}


///////////////////////////////////////////////////////////
//////////////////// Reference Pose ///////////////////////

void
attributes_for_saved_reference_pose(
	utility::tag::AttributeList & attributes,
	std::string const & attribute_name)
{
	attributes_for_saved_reference_pose_w_description( attributes, "", attribute_name );
}

void
attributes_for_saved_reference_pose_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description,
	std::string const & attribute_name)
{
	attributes + utility::tag::XMLSchemaAttribute( attribute_name, utility::tag::xs_string,
		( description == "" ? "Name of reference pose to use" : description ) );
}

core::pose::PoseOP
saved_reference_pose( utility::tag::TagCOP const in_tag, basic::datacache::DataMap & data_map, std::string const & tag_name ){

	if ( in_tag->hasOption(tag_name) ) {
		core::pose::PoseOP refpose(nullptr);
		std::string refpose_name(in_tag->getOption<std::string>( tag_name) );
		TR<<"Loading PDB: "<<refpose_name<<std::endl;

		if ( !data_map.has("spm_ref_poses",refpose_name) ) {
			refpose = core::pose::PoseOP( new core::pose::Pose() );
			data_map.add("spm_ref_poses",refpose_name,refpose );
		} else refpose = data_map.get_ptr<core::pose::Pose>("spm_ref_poses",refpose_name );

		return refpose;
	} else std::cerr << "WARNING: saved_reference_pose function called even though tag has no " + tag_name + " entry. something's unclean somewhere." << std::endl;
	return nullptr;
}


/////////////////////////////////////////////////////////
//////////////////// MoveMap ////////////////////////////

/// @brief utility function for parse_movemap which goes over each of the tags in a movemap section
void
foreach_movemap_tag(
	utility::tag::TagCOP const in_tag,
	core::select::movemap::MoveMapFactoryOP mmf
){
	using namespace utility::tag;
	using namespace core::select;
	using namespace core::select::movemap;

	for ( TagCOP const tag : in_tag->getTags() ) {
		std::string const name( tag->getName() );
		runtime_assert( name == "Jump" || name == "Chain" || name == "Span" );
		if ( name == "Jump" ) {
			auto const num( tag->getOption< int >( "number" ) );
			bool const setting( tag->getOption< bool >( "setting" ) );
			if ( num == 0 ) { // set all jumps if number==0
				mmf->all_jumps( setting );
			} else {
				move_map_action action( setting ? mm_enable : mm_disable );
				// We set soft_error here, as existing usage can sometime set jumps that don't exist in the pose.
				jump_selector::JumpSelectorCOP jump( new jump_selector::JumpIndexSelector( num, /* soft_error= */ true ) );
				mmf->add_jump_action( action, jump );
			}
		}
		if ( name == "Chain" ) {
			std::string const chain( tag->getOption< std::string >( "number" ) );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			bool const bondangle( tag->getOption< bool >( "bondangle", false ) );
			bool const bondlength( tag->getOption< bool >( "bondlength", false ) );

			residue_selector::ResidueSelectorCOP chain_select( new residue_selector::ChainSelector( chain ) );

			mmf->add_bb_action( (bb ? mm_enable : mm_disable), chain_select );
			mmf->add_chi_action( (chi ? mm_enable : mm_disable), chain_select );
			if ( bondangle || bondlength ) {
				mmf->add_bondangles_action( (bondlength ? mm_enable : mm_disable), chain_select );
				mmf->add_bondlengths_action( (bondlength ? mm_enable : mm_disable), chain_select );
			}
		}
		if ( name == "Span" ) {
			std::string const begin( tag->getOption< std::string >( "begin" ) );
			std::string const end( tag->getOption< std::string >( "end" ) );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			bool const bondangle( tag->getOption< bool >( "bondangle", false ) );
			bool const bondlength( tag->getOption< bool >( "bondlength", false ) );

			runtime_assert( ! begin.empty() );
			runtime_assert( ! end.empty() );
			// For backward compatibility, ResidueSpan selector allows ranges outside the pose
			residue_selector::ResidueSelectorCOP span_select( new residue_selector::ResidueSpanSelector( begin, end ) );

			mmf->add_bb_action( (bb ? mm_enable : mm_disable), span_select );
			mmf->add_chi_action( (chi ? mm_enable : mm_disable), span_select );
			if ( bondangle || bondlength ) {
				mmf->add_bondangles_action( (bondlength ? mm_enable : mm_disable), span_select );
				mmf->add_bondlengths_action( (bondlength ? mm_enable : mm_disable), span_select );
			}
		}
	}//for tag
}

void
parse_movemap_tag(
	utility::tag::TagCOP const in_tag,
	core::select::movemap::MoveMapFactoryOP mmf
) {
	if ( in_tag->hasOption( "bb" ) ) mmf->all_bb( in_tag->getOption< bool >( "bb" ) );
	if ( in_tag->hasOption( "chi" ) ) mmf->all_chi( in_tag->getOption< bool >( "chi" ) );
	if ( in_tag->hasOption( "jump" ) ) mmf->all_jumps( in_tag->getOption< bool >( "jump" ) );
}


core::select::movemap::MoveMapFactoryOP
parse_movemap_factory_legacy(
	utility::tag::TagCOP in_tag,
	basic::datacache::DataMap & data,
	bool const reset_movemap, /*= true // should we turn everything to true at start?*/
	core::select::movemap::MoveMapFactoryOP mmf_to_modify /* = nullptr */)
{
	using utility::tag::TagCOP;

	TagCOP movemap_tag;

	if ( in_tag->getName() != "MoveMap" ) {
		utility::vector1< TagCOP > const branch_tags( in_tag->getTags() );
		utility::vector1< TagCOP >::const_iterator tag_it;
		for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
			if ( (*tag_it)->getName() =="MoveMap" ) {
				break;
			}
		}
		if ( tag_it == branch_tags.end() ) return mmf_to_modify; // No MoveMap specification in tag
		movemap_tag = *tag_it;
	} else {
		movemap_tag = in_tag;
	}

	core::select::movemap::MoveMapFactoryOP mmf( mmf_to_modify );
	bool from_datamap( false );

	// If we already have the factory, then make a copy for potential further alteration.
	// (I'd prefer to return it directly, but the legacy code allowed further alterations.
	// The one thing I'm *not* doing is re-injecting the alterations back into the datamap,
	// which is something the legacy code (inadvertantly?) did.)
	if ( movemap_tag->hasOption("name") ) {
		std::string const name( movemap_tag->getOption< std::string >( "name" ) );
		if ( data.has( core::select::movemap::movemap_factory_category(), name ) ) {
			from_datamap = true;
			if ( mmf ) {
				TR.Warning << "Ignoring 'default' MoveMapFactory in parse_movemap_factory_legacy(), using existing MoveMapFactory " << name << " instead" << std::endl;
			}
			mmf = core::select::movemap::get_movemap_factory(name, data)->clone(); // Clone!
		}
	}
	if ( mmf == nullptr ) {
		mmf = core::select::movemap::MoveMapFactoryOP( new core::select::movemap::MoveMapFactory );
	}

	if ( reset_movemap && ! from_datamap ) {
		mmf->all_bb( true );
		mmf->all_chi( true );
		mmf->all_jumps( true );
	}

	parse_movemap_tag( movemap_tag, mmf );
	foreach_movemap_tag( movemap_tag, mmf );

	if ( ! from_datamap && movemap_tag->hasOption("name") ) {
		std::string const name( movemap_tag->getOption< std::string >( "name" ) );
		data.add( core::select::movemap::movemap_factory_category(), name, mmf->clone() ); // Store clone, such that return value can be modified
	}

	return mmf;
}

core::select::movemap::MoveMapFactoryOP
parse_named_movemap_factory_legacy(
	utility::tag::TagCOP in_tag,
	std::string const & name,
	basic::datacache::DataMap & data,
	bool const reset_movemap, /*= true // should we turn everything to true at start?*/
	core::select::movemap::MoveMapFactoryOP mmf_to_modify /* = nullptr */)
{
	using utility::tag::TagCOP;

	utility::vector1< TagCOP > const branch_tags( in_tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() =="MoveMap" &&
				(*tag_it)->hasOption("name") &&
				(*tag_it)->getOption< std::string >( "name" ) == name ) {
			break;
		}
	}
	if ( tag_it == branch_tags.end() ) return nullptr; // No MoveMapFactory in tag

	return parse_movemap_factory_legacy(*tag_it, data, reset_movemap, mmf_to_modify);

}

std::string move_map_tag_namer( std::string const & subelement_name ) { return "move_map_" + subelement_name + "_type"; }
std::string optionally_named_move_map_ct_namer( std::string const & ) { return "optionally_named_move_map_type"; }

void
common_movemap_complex_type_def( utility::tag::XMLSchemaComplexTypeGenerator & ct_gen )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList movemap_subelements;
	AttributeList jump_attributes;
	jump_attributes
		+ XMLSchemaAttribute::required_attribute( "number",  xsct_non_negative_integer , "Which jump number (in the FoldTree)" )
		+ XMLSchemaAttribute::required_attribute( "setting", xsct_rosetta_bool , "true for move, false for don't move" );
	movemap_subelements.add_simple_subelement( "Jump", jump_attributes , "jumps are the not-chemistry internal coordinate connections between separate parts of your pose");

	std::string const chi_desc("move sidechain chi torsions?");
	std::string const bb_desc("move backbone torsions?");

	AttributeList chain_attributes;
	chain_attributes
		+ XMLSchemaAttribute::required_attribute( "number", xsct_non_negative_integer , "which chain?" )
		+ XMLSchemaAttribute::required_attribute( "chi",    xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute::required_attribute( "bb",     xsct_rosetta_bool , bb_desc );
	movemap_subelements.add_simple_subelement( "Chain", chain_attributes , "this controls a kinematically contiguous chain (think protein chains)");

	AttributeList span_attributes;
	span_attributes
		+ XMLSchemaAttribute::required_attribute( "begin",     xsct_non_negative_integer , "beginning of span" )
		+ XMLSchemaAttribute::required_attribute( "end",       xsct_non_negative_integer , "end of span" )
		+ XMLSchemaAttribute::required_attribute( "chi",       xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute::required_attribute( "bb",        xsct_rosetta_bool , bb_desc )
		+ XMLSchemaAttribute( "bondangle", xsct_rosetta_bool , "move 3-body angles?" )
		+ XMLSchemaAttribute( "bondlength", xsct_rosetta_bool , "move 2-body lengths?" );
	movemap_subelements.add_simple_subelement( "Span", span_attributes , "XRW TO DO, probably a user-defined region of the Pose");
	movemap_subelements.complex_type_naming_func( & move_map_tag_namer );


	AttributeList movemap_tag_attributes;
	movemap_tag_attributes
		+ XMLSchemaAttribute( "bb", xsct_rosetta_bool , bb_desc )
		+ XMLSchemaAttribute( "chi", xsct_rosetta_bool , chi_desc )
		+ XMLSchemaAttribute( "jump", xsct_rosetta_bool , "move all jumps?" );

	ct_gen.element_name( "MoveMap" )
		.description( "MoveMaps dicate what degrees of freedom are mobile to other bits of code, especially minimization - they do NOT work with packing" )
		.add_attributes( movemap_tag_attributes )
		.set_subelements_repeatable( movemap_subelements );

}

void
append_subelement_for_parse_movemap_factory_legacy(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::XMLSchemaSimpleSubelementList & subelements,
	std::string const & description
)
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGenerator named_move_map;
	common_movemap_complex_type_def( named_move_map );
	named_move_map
		.description( description == "" ? "MoveMap specification" : description )
		.add_optional_name_attribute()
		.complex_type_naming_func( & optionally_named_move_map_ct_namer )
		.write_complex_type_to_schema( xsd );

	subelements.add_already_defined_subelement( "MoveMap", & optionally_named_move_map_ct_namer );
}

/////////////////////////////////////////////////////////
//////////////////// Filter ////////////////////////////

protocols::filters::FilterOP
parse_filter( std::string const & filter_name, protocols::filters::Filters_map const & filters ){
	auto filter_it( filters.find( filter_name ) );
	if ( filter_it == filters.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Filter "+filter_name+" not found" );
	}
	return filter_it->second;
}



/////////////////////////////////////////////////////////
//////////////////// Mover //////////////////////////////

protocols::moves::MoverOP
parse_mover( std::string const & mover_name, protocols::moves::Movers_map const & movers ){
	auto mover_it( movers.find( mover_name ) );
	if ( mover_it == movers.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Mover "+mover_name+" not found" );
	}
	return mover_it->second;
}



/////////////////////////////////////////////////////////
//////////////////// XYZ Vector /////////////////////////

void
attributes_for_parse_xyz_vector( utility::tag::AttributeList & attlist ){
	using namespace utility::tag;
	attlist
		+ XMLSchemaAttribute::required_attribute( "x", xsct_real, "X coordinate for this vector" )
		+ XMLSchemaAttribute::required_attribute( "y", xsct_real, "Y coordinate for this vector" )
		+ XMLSchemaAttribute::required_attribute( "z", xsct_real, "Z coordinate for this vector" );
}


/// @brief utility function for parsing xyzVector
numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagCOP const xyz_vector_tag ){
	if ( ! xyz_vector_tag->hasOption("x") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "xyz_vector requires 'x' coordinates option");
	if ( ! xyz_vector_tag->hasOption("y") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "xyz_vector requires 'y' coordinates option");
	if ( ! xyz_vector_tag->hasOption("z") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "xyz_vector requires 'z' coordinates option");

	numeric::xyzVector< core::Real > xyz_v (
		xyz_vector_tag->getOption<core::Real>("x"),
		xyz_vector_tag->getOption<core::Real>("y"),
		xyz_vector_tag->getOption<core::Real>("z")
	);

	return xyz_v;

}


// void
// attributes_for_report_to_db( utility::tag::AttributeList & attlist, utility::tag::XMLSchemaDefinition & xsd)
// {
//  using namespace utility::tag;
//
//  XMLSchemaRestriction database_mode;
//  database_mode.name( "database_mode_string" );
//  database_mode.base_type( xs_string );
//  database_mode.add_restriction( xsr_enumeration, "sqlite3" );
//  database_mode.add_restriction( xsr_enumeration, "postgres" );
//  database_mode.add_restriction( xsr_enumeration, "mysql" );
//  xsd.add_top_level_element( database_mode );
//  XMLSchemaRestriction database_transaction_mode;
//  database_transaction_mode.name( "database_transaction_mode_string" );
//  database_transaction_mode.base_type( xs_string );
//  database_transaction_mode.add_restriction( xsr_enumeration, "none" );
//  database_transaction_mode.add_restriction( xsr_enumeration, "standard" );
//  database_transaction_mode.add_restriction( xsr_enumeration, "chunk" );
//  xsd.add_top_level_element( database_transaction_mode );
//
//
//
//  XMLSchemaRestriction mode;
//  mode.name( "report_to_db_relevant_residues_mode" );
//  mode.base_type( xs_string );
//  //can be "explicit" or "implicit" with any capitalization
//  mode.add_restriction( xsr_pattern, "[eEiI][xXmM][pP][lL][iI][cC][iI][tT]" );
//  xsd.add_top_level_element( mode );
//
//
//  attlist
//   + XMLSchemaAttribute( "resource_description", xs_string, "Resource description referring to the resource to be output" )
//   + XMLSchemaAttribute( "database_resource", xs_string, "Resource description referring to the output database" )
//   + XMLSchemaAttribute( "database_resource_tag", xs_string, "Tag referring to the output database" )
//   + XMLSchemaAttribute::attribute_w_default( "transaction_mode", "database_transaction_mode_string", "Transaction mode for database output", "standard" )
//   + XMLSchemaAttribute( "database_mode", "database_mode_string", "Which type of output database to use?" )
//   + XMLSchemaAttribute( "database_name", xs_string, "Name of output database" )
//   + XMLSchemaAttribute( "database_pq_schema", xs_string, "Schema name within the database" )
//   + XMLSchemaAttribute( "database_separate_db_per_mpi_process", xsct_rosetta_bool, "Use a separate database for each MPI process? Specific to sqlite3" )
//   + XMLSchemaAttribute( "database_partition", xs_integer, "Database partition to use" )
//   + XMLSchemaAttribute( "database_read_only", xsct_rosetta_bool, "Is the database read-only?" )
//   + XMLSchemaAttribute( "database_host", xs_string, "URI to the database server (postgres only)" )
//   + XMLSchemaAttribute( "database_user", xs_string, "Username for database access( postgres only)" )
//   + XMLSchemaAttribute( "database_password", xs_string, "Password for database access (postgres only)" )
//   + XMLSchemaAttribute( "database_port", xsct_non_negative_integer, "Port to use for access to database server (postgres only)" )
//   + XMLSchemaAttribute( "batch_description", xs_string, "Description of features database" )
//   + XMLSchemaAttribute( "protocol_id", xsct_non_negative_integer, "Manually controls protocol_id associated with this ReportToDB tag. Autoincrements by default." )
//   + XMLSchemaAttribute( "batch_id", xsct_non_negative_integer, "Manually controls the batch_id associated with this ReportToDB tag. Autoincrements by default." )
//   + XMLSchemaAttribute( "use_transactions", xsct_rosetta_bool, "Use transactions to group database i/o to be more efficient. Turning them off can help debugging." )
//   + XMLSchemaAttribute::attribute_w_default( "cache_size", xsct_non_negative_integer, "Specify the maximum number 1k pages to keep in memory before writing to disk.", "2000" )
//   + XMLSchemaAttribute::attribute_w_default( "remove_xray_virt", xsct_rosetta_bool, "Remove virtual residue attached during xray refine process", "false" )
//   + XMLSchemaAttribute::attribute_w_default( "relevant_residues_mode", "report_to_db_relevant_residues_mode", "Determine what features are reported given the relevant residues", "explicit" );
//
//  rosetta_scripts::attributes_for_parse_task_operations( attlist );
// }

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
		protocols::rosetta_scripts::RosettaScriptsParser parser = protocols::rosetta_scripts::RosettaScriptsParser();
		schema = parser.xsd_for_rosetta_scripts();
	} catch ( utility::excn::Exception const & e ) {
		TR.Error << "An error was encountered while attempting to generate the XSD schema - it will not be saved.\n";
		TR.Error << e.msg() << "\n";
		return;
	}

	std::ofstream ofs( filename );
	ofs << schema << std::endl;
}


} //RosettaScripts
} //protocols
