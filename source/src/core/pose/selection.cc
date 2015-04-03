// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/selection.cc
/// @brief pose residue selections
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Jacob Corn (jecorn@u.washington.edu)
///	@author Rocco Moretti (rmoretti@u.washington.edu)
/// @author Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <core/pose/selection.hh>
// Project Headers

#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <boost/foreach.hpp>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/string.functions.hh>


static thread_local basic::Tracer TR( "core.pose.selection" );

namespace core {
namespace pose {

using namespace core::scoring;
using namespace core;
using namespace std;
using utility::vector1;

/// @brief a convenience function to test whether the user has specified pdb numbering rather than rosetta numbering.
core::Size
get_resnum( utility::tag::TagCOP tag_ptr, core::pose::Pose const & pose, std::string const & prefix/*=""*/ ) {
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
parse_resnum(
	std::string const& resnum,
	core::pose::Pose const& pose)
{

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
		runtime_assert(info != 0);
		return info->pdb2pose( chain[0], n );
	}
	else { // Rosetta Number
		TR.Trace << "Interpreting " << n << " as a Rosetta residue number." << std::endl;
		return n;
	}
}

/// @brief Extracts residue numbers from a 'selection'.
/// @detail Recognizes two forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
///   - name3=ALA selects all alanines
/// @return the rosetta residue numbers for the string, or 0 upon an error
utility::vector1<core::Size>
parse_selection_block(
	std::string const& sele,
	core::pose::Pose const& pose
){
	utility::vector1<Size> res;
	if(sele.size()==9&&sele.substr(0,6)=="name3="){
		string name3 = sele.substr(6,3);
		for(Size ir=1; ir <= pose.n_residue(); ++ir) {
			if(pose.residue(ir).name3()==name3){
				res.push_back(ir);
			}
		}
	} else {
		res.push_back(parse_resnum(sele,pose));
	}
	return res;
}


/// @brief Extracts a list of residue numbers from a tag.
/// @details The tag should contain a comma-separated list of numbers, in either
///   pdb or rosetta format (@see parse_resnum for details)
vector1<Size>
get_resnum_list(
	utility::tag::TagCOP tag_ptr,
	string const& tag,
	pose::Pose const& pose
){
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
get_resnum_list(
	std::string const str,
	core::pose::Pose const & pose
){
	using namespace std;
	using namespace utility;
	set< Size > resid;

	resid.clear();
	vector1< string> const str_residues( utility::string_split( str , ',' ) );
	BOOST_FOREACH( string const res, str_residues ){
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
		utility::vector1<core::Size> const nums( parse_selection_block( res, pose ) );
		for(utility::vector1<core::Size>::const_iterator i = nums.begin(); i != nums.end(); ++i){
			Size num = *i;
			runtime_assert( num );
			resid.insert( num );
		}
	}//foreach
	return resid;
}


// fpd same as 'get_resnum_list', but preserve ordering from input list
utility::vector1<Size>
get_resnum_list_ordered(
	std::string const str,
	core::pose::Pose const & pose
){
	using namespace std;
	using namespace utility;
	utility::vector1<Size> resid;

	resid.clear();
	vector1< string> const str_residues( utility::string_split( str , ',' ) );
	BOOST_FOREACH( string const res, str_residues ){
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
					resid.push_back( i );
				continue;
			}
		}
		utility::vector1<core::Size> const nums( parse_selection_block( res, pose ) );
		for(utility::vector1<core::Size>::const_iterator i = nums.begin(); i != nums.end(); ++i){
			Size num = *i;
			runtime_assert( num );
			resid.push_back( num );
		}
	}//foreach
	return resid;
}


} //pose
} //core
