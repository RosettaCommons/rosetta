// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/RosettaScripts/util.cc
/// @brief Utility functions useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu)

// Unit Headers
#include <protocols/rosetta_scripts/util.hh>

// Project Headers

#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// C++ headers

static basic::Tracer TR( "protocols.RosettaScripts.util" );

namespace protocols {
namespace rosetta_scripts {

using namespace core::scoring;
using namespace protocols::moves;
using namespace core;
using namespace std;
using utility::vector1;

// a convenience function to test whether the user has specified pdb numbering rather than rosetta numbering.
core::Size
get_resnum( utility::tag::TagPtr const tag_ptr, core::pose::Pose const & pose, std::string const & prefix/*=""*/ ) {
	core::Size resnum( 0 );
	bool const pdb_num_used( tag_ptr->hasOption( prefix + "pdb_num" ) );
	if( pose.pdb_info().get() == NULL ){//no pdbinfo for this pose (e.g., silent file), resort to using just the residue number
		if( pdb_num_used ){
			TR<<"Bad tag: "<< *tag_ptr<<std::endl;
			utility_exit_with_message( "pdb_num used but no pdb_info found. Use res_num instead" );
			return( 0 );
		}
	}
	else{
		core::pose::PDBPoseMap const pose_map( pose.pdb_info()->pdb2pose() );
		if( pdb_num_used ) {
			std::string pdbnum( tag_ptr->getOption<std::string>( prefix + "pdb_num" ) );
			char const chain( pdbnum[ pdbnum.length() - 1 ] );
			std::stringstream ss( pdbnum.substr( 0, pdbnum.length() - 1 ) );
			core::Size number;
			ss >> number;
			resnum = pose_map.find( chain, number );
		}
	}
	if( !pdb_num_used )
		resnum = tag_ptr->getOption<core::Size>( prefix + "res_num" );

	runtime_assert( resnum );
	return( resnum );
}

/// @brief Extracts a residue number from a string.
/// @detail Recognizes two forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
/// @return the rosetta residue number for the string, or 0 upon an error
core::Size
parse_resnum(std::string const& resnum, core::pose::Pose const& pose) {

	string::const_iterator input_end = resnum.end();
	//Set number to the sequence of digits at the start of input [0-9]*
	string::const_iterator number_start = resnum.begin();
	string::const_iterator number_end = resnum.begin();
	while( number_end != input_end && *number_end >= '0' && *number_end <= '9' ) {
		++number_end;
	}
	//Set chain to the following characters
	string::const_iterator chain_start = number_end;
	string::const_iterator chain_end = number_end;
	while(  chain_end != input_end
		&& (('A' <= *chain_end && *chain_end <= 'Z') ||
			('a' <= *chain_end && *chain_end <= 'z') ||
			'_' == *chain_end ) )
	{
		++chain_end;
	}

	string number(number_start,number_end);
	string chain(chain_start,chain_end);

	//Require that the whole string match, and that the chain be a single char
	if( chain_end != input_end || chain.size() > 1 || number.size() < 1) {
		TR.Error << "Could not parse '" << resnum << "' into a residue number." << std::endl;
		return Size(0);
	}

	Size n;
	std::istringstream ss( number );
	ss >> n;
	if( chain.size() == 1 ) { // PDB Number
		TR.Trace << "Interpretting " << n << chain << " as a pdb number." << std::endl;
		pose::PDBInfoCOP info = pose.pdb_info();
		runtime_assert(info);
		return info->pdb2pose( chain[0], n );
	}
	else { // Rosetta Number
		TR.Trace << "Interpreting " << n << " as a Rosetta residue number." << std::endl;
		return n;
	}
}


/// @brief Extracts a list of residue numbers from a tag.
/// @details The tag should contain a comma-separated list of numbers, in either
///   pdb or rosetta format (@see parse_resnum for details)
vector1<Size>
get_resnum_list(utility::tag::TagPtr const tag_ptr, string const& tag, pose::Pose const& pose)
{
	vector1< Size > resnums;
	if( ! tag_ptr->hasOption( tag ) ) {
		TR<<"Error: No "<<tag<<" option was found in tag "<<tag_ptr<<std::endl;
		utility_exit();
		return resnums;
	}

	set<Size> const resnums_set( get_resnum_list( tag_ptr->getOption< std::string >( tag ), pose ) );
	resnums.clear();
	resnums.insert( resnums.begin(), resnums_set.begin(), resnums_set.end() );
	sort( resnums.begin(), resnums.end() );
	unique( resnums.begin(), resnums.end() );

	return resnums;
}

set<Size>
get_resnum_list( std::string const str, core::pose::Pose const & pose )
{
	using namespace std;
	using namespace utility;
	set< Size > resid;

	resid.clear();
	vector1< string> const str_residues( utility::string_split( str , ',' ) );
	foreach( string const res, str_residues ){
		if( res == "" ) continue;
		if( res.find('-') != string::npos) {
			// Handle residue range
			vector1< string> const str_limits( utility::string_split( res , '-' ) );
			if ( str_limits.size() == 2) {
				core::Size const start ( parse_resnum( str_limits[1], pose ) );
				core::Size const end ( parse_resnum( str_limits[2], pose ) );
				if ( start && end && start > end ) {
					utility_exit_with_message("Invalid residue range: " + res);
				}
				for(core::Size i = start; i <= end; ++i )
					resid.insert( i );
				continue;
			}
		}
		core::Size const num( parse_resnum( res, pose ) );
		runtime_assert( num );
		resid.insert( num );
	}//foreach
	return resid;
}


/// @details This is essentially a shameless copy of Justin's PackRotamersMover::parse_task_operations. In truth
/// DesignRepackMover should disappear into Justin's better organized class, but this will wait... (SJF)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

  TaskFactoryOP new_task_factory( new TaskFactory );
  if ( ! tag->hasOption("task_operations") ) return new_task_factory;
  std::string const t_o_val( tag->getOption<std::string>("task_operations") );
  typedef utility::vector1< std::string > StringVec;
  StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations to mover "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";
  for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
        t_o_key != end; ++t_o_key ) {
    if ( data.has( "task_operations", *t_o_key ) ) {
      new_task_factory->push_back( data.get< TaskOperation * >( "task_operations", *t_o_key ) );
			TR<<*t_o_key<<' ';
    } else {
      utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
    }
  }
	TR<<std::endl;
  return new_task_factory;
}

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data )
{
	using core::pack::task::operation::TaskOperationOP;
	using core::pack::task::operation::TaskOperation;
	typedef std::string String;

	utility::vector1< TaskOperationOP > task_operations;
	String const t_o_val( tag->getOption<String>("task_operations", "" ) );
	if( t_o_val != "" ){
		utility::vector1< String > const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( utility::vector1< String >::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
				t_o_key != end; ++t_o_key ) {
			if ( data.has( "task_operations", *t_o_key ) ) {
				task_operations.push_back( data.get< TaskOperation* >( "task_operations", *t_o_key ) );
			} else {
				utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
			}
		}
	}
	return task_operations;
}

/// @details Utility function to find a scorefunction from parser-provided data. This is essentially a shameless
/// copy of Justin's PackRotamersMover::parse_score_function.
core::scoring::ScoreFunctionOP
parse_score_function( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data )
{
	std::string const scorefxn_key( tag->getOption<std::string>("scorefxn", "score12" ) );
	if ( ! data.has( "scorefxns", scorefxn_key ) ) {
		utility_exit_with_message("ScoreFunction " + scorefxn_key + " not found in DataMap.");
	}
	return data.get< ScoreFunction* >( "scorefxns", scorefxn_key );
}

core::pose::PoseOP
saved_reference_pose( utility::tag::TagPtr const in_tag, protocols::moves::DataMap & data_map ){

	if( in_tag->hasOption("reference_name") ){
		core::pose::PoseOP refpose(NULL);
		std::string refpose_name(in_tag->getOption<std::string>( "reference_name") );

		if( !data_map.has("spm_ref_poses",refpose_name) ){
			refpose = new core::pose::Pose();
			data_map.add("spm_ref_poses",refpose_name,refpose );
		}
		else refpose = data_map.get<core::pose::Pose *>("spm_ref_poses",refpose_name );

		return refpose;
	}
	else std::cerr << "WARNING: saved_reference_pose function called even though tag has no 'reference_name' entry. something's unclean somewhere." << std::endl;
	return NULL;
}

/// @brief utility function for parse_movemap which goes over each of the tags in a movemap section
void
foreach_movemap_tag( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm ){
	using namespace core::kinematics;
	using namespace utility::tag;

	foreach( TagPtr const tag, in_tag->getTags() ){
		std::string const name( tag->getName() );
		runtime_assert( name == "Jump" || name == "Chain" || name == "Span" );
		if( name == "Jump" ){
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const setting( tag->getOption< bool >( "setting" ) );
			if( num == 0 ) mm->set_jump( setting ); // set all jumps if number==0
			else mm->set_jump( num, setting );
		}
		if( name == "Chain" ){
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			core::Size const chain_begin( pose.conformation().chain_begin( num ) );
			core::Size const chain_end( pose.conformation().chain_end( num ) );
			for( core::Size i( chain_begin ); i <= chain_end; ++i ){
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
		}
		if( name == "Span" ){
			core::Size const begin( tag->getOption< core::Size >( "begin" ) );
			core::Size const end( tag->getOption< core::Size >( "end" ) );
			runtime_assert( end >= begin );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			for( core::Size i( begin ); i <= end; ++i ){
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
		}
	}//foreach tag
}

/// @brief variant of parse_movemap that takes in a datamap and searches it for already existing movemaps
void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm, protocols::moves::DataMap & data, bool const reset_map ){
	using utility::tag::TagPtr;
	using namespace core::kinematics;

	if( in_tag() == NULL ) return;

	utility::vector1< TagPtr > const branch_tags( in_tag->getTags() );
	utility::vector1< TagPtr >::const_iterator tag_it;
	for( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ){
		if( (*tag_it)->getName() == "MoveMap" ){
			break;
		}
	}
	if( reset_map ){
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}
	if( tag_it == branch_tags.end() ) return;

	if( (*tag_it)->hasOption( "name" ) ){
		std::string const name( (*tag_it)->getOption< std::string >( "name" ) );
		if( data.has( "movemaps", name ) ){
			mm = data.get< MoveMap * >( "movemaps", name );
			TR<<"Found movemap "<<name<<" on datamap"<<std::endl;
		}
		else{
			data.add( "movemaps", name, mm );
			TR<<"Adding movemap "<<name<<" to datamap"<<std::endl;
		}
	}
	foreach_movemap_tag( *tag_it, pose, mm );
}

///@details modifies an existing movemap according to tag
/// the movemap defaults to move all bb, chi, and jumps.
void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm ){
	using utility::tag::TagPtr;
	using namespace core::kinematics;

	if( in_tag() == NULL ) return;

	utility::vector1< TagPtr > const branch_tags( in_tag->getTags() );
	utility::vector1< TagPtr >::const_iterator tag_it;
	for( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ){
		if( (*tag_it)->getName() == "MoveMap" ){
			break;
		}
	}
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( true );
	if( tag_it == branch_tags.end() ) return;

	if( (*tag_it)->hasOption( "name" ) ){
		TR<<"ERROR in "<<*tag_it<<'\n';
		utility_exit_with_message( "Tag called with option name but this option is not available to this mover. Note that FastRelax cannot work with a prespecified movemap b/c its movemap is parsed at apply time. Sorry." );
	}
	if( tag_it == branch_tags.end() ) return;

	foreach_movemap_tag( *tag_it, pose, mm );
}

protocols::filters::FilterOP
parse_filter( std::string const filter_name, protocols::filters::Filters_map const & filters ){
  protocols::filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );
  if( filter_it == filters.end() )
    utility_exit_with_message( "Filter "+filter_name+" not found" );
  return filter_it->second;
}

protocols::moves::MoverOP
parse_mover( std::string const mover_name, protocols::moves::Movers_map const & movers ){
  protocols::moves::Movers_map::const_iterator mover_it( movers.find( mover_name ) );
  if( mover_it == movers.end() )
    utility_exit_with_message( "Mover "+mover_name+" not found" );
  return mover_it->second;
}

/// @brief utility function for parsing xyzVector
numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagPtr const xyz_vector_tag ){
	if ( ! xyz_vector_tag->hasOption("x") ) utility_exit_with_message("xyz_vector requires 'x' coordinates option");
	if ( ! xyz_vector_tag->hasOption("y") ) utility_exit_with_message("xyz_vector requires 'y' coordinates option");
	if ( ! xyz_vector_tag->hasOption("z") ) utility_exit_with_message("xyz_vector requires 'z' coordinates option");

	numeric::xyzVector< core::Real > xyz_v (
		xyz_vector_tag->getOption<core::Real>("x"),
		xyz_vector_tag->getOption<core::Real>("y"),
		xyz_vector_tag->getOption<core::Real>("z")
	);

	return xyz_v;

}

} //RosettaScripts
} //protocols
