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
/// @author Rocco Moretti (rmoretti@u.washington.edu)
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


static THREAD_LOCAL basic::Tracer TR( "core.pose.selection" );

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
	if ( pose.pdb_info().get() == NULL ) { //no pdbinfo for this pose (e.g., silent file), resort to using just the residue number
		if ( pdb_num_used ) {
			TR<<"Bad tag: "<< *tag_ptr<<std::endl;
			utility_exit_with_message( "pdb_num used but no pdb_info found. Use res_num instead" );
			return( 0 );
		}
	} else {
		core::pose::PDBPoseMap const pose_map( pose.pdb_info()->pdb2pose() );
		if ( pdb_num_used ) {
			std::string pdbnum( tag_ptr->getOption<std::string>( prefix + "pdb_num" ) );
			char const chain( pdbnum[ pdbnum.length() - 1 ] );
			std::stringstream ss( pdbnum.substr( 0, pdbnum.length() - 1 ) );
			core::Size number(0);
			ss >> number;
			resnum = pose_map.find( chain, number );
		}
	}
	if ( !pdb_num_used ) {
		resnum = tag_ptr->getOption<core::Size>( prefix + "res_num" );
	}

	runtime_assert( resnum );
	return( resnum );
}

/// @brief Extracts a residue number from a string.
/// @detail Recognizes three forms of numbering:
///   - Rosetta residue numbers (numbered sequentially from 1 to the last residue
///     in the pose). These have the form [0-9]+
///   - PDB numbers. These have the form [0-9]+[A-Z], where the trailing letter
///     is the chain ID.
///  - Reference pose numbers.  These have the form refpose([refpose name],[refpose number]).
/// In addition, relative numbers are permitted (of the form +[number] or -[number]) in conjunction
/// with reference pose numbering.  For example, one might say "refpose(state1,17)+3", which means
/// three residues past the residue correpsonding to residue 17 in the reference pose called "state1".
/// @return the rosetta residue number for the string, or 0 upon an error
core::Size
parse_resnum(
	std::string const& resnum,
	core::pose::Pose const& pose,
	bool const check_for_refpose
) {
	// Scope 1: Check whether this string is of the form
	// "resnum(<refpose_name>,<refpose_number>)+/-<offset>"
	// and parse accordingly if it is.
	if ( check_for_refpose ) {
		std::string refpose_name("");
		core::Size refpose_number(0);
		signed long refpose_offset(0); //Must be signed!
		if ( is_referencepose_number(resnum, refpose_name, refpose_number, refpose_offset) ) {
			return get_resnumber_from_reference_pose(refpose_name, refpose_number, refpose_offset, pose);
		}
	}
	// Otherwise, this is NOT of that form, and we should parse it as a Rosetta
	// number or as a PDB number:

	string::const_iterator input_end = resnum.end();
	//Set number to the sequence of digits at the start of input [0-9]*
	string::const_iterator number_start = resnum.begin();
	string::const_iterator number_end = resnum.begin();
	while ( number_end != input_end && *number_end >= '0' && *number_end <= '9' ) {
		++number_end;
	}
	//Set chain to the following characters
	string::const_iterator chain_start = number_end;
	string::const_iterator chain_end = number_end;
	while (  chain_end != input_end
			&& (('A' <= *chain_end && *chain_end <= 'Z') ||
			('a' <= *chain_end && *chain_end <= 'z') ||
			'_' == *chain_end ) )
			{
		++chain_end;
	}

	string number(number_start,number_end);
	string chain(chain_start,chain_end);

	//Require that the whole string match, and that the chain be a single char
	if ( chain_end != input_end || chain.size() > 1 || number.size() < 1 ) {
		TR.Error << "Could not parse '" << resnum << "' into a residue number." << std::endl;
		return Size(0);
	}

	Size n(0);
	std::istringstream ss( number );
	ss >> n;
	if ( chain.size() == 1 ) { // PDB Number
		TR.Trace << "Interpretting " << n << chain << " as a pdb number." << std::endl;
		pose::PDBInfoCOP info = pose.pdb_info();
		runtime_assert(info != 0);
		return info->pdb2pose( chain[0], n );
	} else { // Rosetta Number
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
	if ( sele.size()==9&&sele.substr(0,6)=="name3=" ) {
		string name3 = sele.substr(6,3);
		for ( Size ir=1; ir <= pose.n_residue(); ++ir ) {
			if ( pose.residue(ir).name3()==name3 ) {
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
	if ( ! tag_ptr->hasOption( tag ) ) {
		TR<<"Error: No "<<tag<<" option was found in tag "<<tag_ptr<<std::endl;
		utility_exit();
		return resnums;
	}
	set<Size> const resnums_set( get_resnum_list( tag_ptr->getOption< std::string >( tag ), pose ) );
	resnums.clear();
	resnums.insert( resnums.begin(), resnums_set.begin(), resnums_set.end() );
	sort( resnums.begin(), resnums.end() );
	vector1<Size>::iterator last = unique( resnums.begin(), resnums.end() );
	resnums.erase(last, resnums.end() );

	return resnums;
}

set<Size>
get_resnum_list(
	std::string const & str,
	core::pose::Pose const & pose
){
	using namespace std;
	using namespace utility;
	set< Size > resid;

	resid.clear();
	vector1< string> const str_residues( utility::string_split( str , ',' ) );
	BOOST_FOREACH ( string const res, str_residues ) {
		if ( res == "" ) continue;
		if ( res.find('-') != string::npos ) {
			// Handle residue range
			vector1< string> const str_limits( utility::string_split( res , '-' ) );
			if ( str_limits.size() == 2 ) {
				core::Size const start ( parse_resnum( str_limits[1], pose ) );
				core::Size const end ( parse_resnum( str_limits[2], pose ) );
				if ( start && end && start > end ) {
					utility_exit_with_message("Invalid residue range: " + res);
				}
				for ( core::Size i = start; i <= end; ++i ) {
					resid.insert( i );
				}
				continue;
			}
		}
		utility::vector1<core::Size> const nums( parse_selection_block( res, pose ) );
		for ( utility::vector1<core::Size>::const_iterator i = nums.begin(); i != nums.end(); ++i ) {
			Size num = *i;
			runtime_assert_string_msg( num, "Error in core::pose::selection::get_resnum_list(): A residue index could not be parsed." );
			resid.insert( num );
		}
	}//foreach
	return resid;
}


// fpd same as 'get_resnum_list', but preserve ordering from input list
utility::vector1<Size>
get_resnum_list_ordered(
	std::string const & str,
	core::pose::Pose const & pose
){
	using namespace std;
	using namespace utility;
	utility::vector1<Size> resid;

	resid.clear();
	vector1< string> const str_residues( utility::string_split( str , ',' ) );
	BOOST_FOREACH ( string const res, str_residues ) {
		if ( res == "" ) continue;
		if ( res.find('-') != string::npos ) {
			// Handle residue range
			vector1< string> const str_limits( utility::string_split( res , '-' ) );
			if ( str_limits.size() == 2 ) {
				core::Size const start ( parse_resnum( str_limits[1], pose ) );
				core::Size const end ( parse_resnum( str_limits[2], pose ) );
				if ( start && end && start > end ) {
					utility_exit_with_message("Invalid residue range: " + res);
				}
				for ( core::Size i = start; i <= end; ++i ) {
					resid.push_back( i );
				}
				continue;
			}
		}
		utility::vector1<core::Size> const nums( parse_selection_block( res, pose ) );
		for ( utility::vector1<core::Size>::const_iterator i = nums.begin(); i != nums.end(); ++i ) {
			Size num = *i;
			runtime_assert( num );
			resid.push_back( num );
		}
	}//foreach
	return resid;
}

/// @brief Is a string of the format "refpose(<refposename>,<refposenumber>)" or "refpose(<refposename>,<refposenumber>)+/-<number>"?
/// @details  If this successfully determines that this is a string of this format, it populates the refpose_string, refpose_resnumber,
/// and refpose_offset variables with the name of the ReferencePose, the number of the residue in the reference pose, and the +/- offset
/// number parsed from this string.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
bool is_referencepose_number(
	std::string const &str,
	std::string &refpose_string,
	core::Size &refpose_resnumber,
	signed long &refpose_offset
)
{

	//TR << "Input string to is_referencepose_number():" << str << std::endl; //DELETE ME

	// 1.  Does the string start with "refpose("?
	if ( str.substr(0,8)!="refpose(" ) return false;

	//TR << "Check 1 passed -- the string begins with \"refpose(\"." << std::endl; //DELETE ME

	// 2.  Does the string have a comma, and then a bracket?
	core::Size commaposition(0);
	core::Size endbracketposition(0);
	for ( core::Size i=8; i<str.length(); ++i ) {
		if ( commaposition==0 ) {
			if ( str[i] == ',' ) commaposition=i;
		} else {
			if ( str[i] == ')' ) {
				endbracketposition=i;
				break;
			}
		}
	}
	if ( commaposition==0 || commaposition<9 || endbracketposition==0 || endbracketposition <= commaposition+1 ) return false;

	//TR << "Check 2 passed -- the string has a comma following the opening of the parentheses, and a close parenthesis, with stuff in between." << std::endl; //DELETE ME

	std::string const refpose_string_out( str.substr( 8, commaposition-8 ) ); //Store what's between the first parenthesis and the comma.  This is the name of the reference pose.
	//TR << "refpose_string_out=" << refpose_string_out << std::endl; //DELETE ME

	//Parse the number in the reference pose:
	std::string const refpose_resnumber_string_out( str.substr(commaposition+1, endbracketposition-commaposition-1) );
	//TR << "refpose_resnumber_string_out=" << refpose_resnumber_string_out << std::endl;
	for ( core::Size i=0, imax=refpose_resnumber_string_out.length(); i<imax; ++i ) {
		if ( refpose_resnumber_string_out[i]<'0' || refpose_resnumber_string_out[i]>'9' ) return false; //Fail if there isn't a number between the comma and the close of the parentheses
	}
	core::Size refpose_resnumber_out(0);
	{ //Scope for temporary istringstream variable:
		std::istringstream ss( refpose_resnumber_string_out );
		ss >> refpose_resnumber_out; //Parse the number into a core::Size.
	}

	//TR << "Check 3 passed -- the number in the reference pose is, indeed, a number, which we have parsed as " << refpose_resnumber_out << "." << std::endl; //DELETE ME

	bool hasoffset(false);
	if ( str.length()>endbracketposition+1 ) { //If the string doesn't end with the end bracket position, then there must be an offset.
		//Check for an offset:
		if ( str[endbracketposition+1]!='+' && str[endbracketposition+1]!='-' ) return false; //Fail if the position after the bracket isn't a + or a -.
		for ( core::Size i=endbracketposition+2, imax=str.length(); i<imax; ++i ) {
			if ( str[i] < '0' || str[i] > '9' ) return false; //Fail if the positions after the +/- contain non-numbers.
		}
		hasoffset = true;
	}

	//TR << "Check 4 passed.  " << (hasoffset ? "There is a properly-formatted offset." : "There is no offset.") << std::endl; //DELETE ME

	//At this point, we know that we will be returning "true" -- the string is formatted properly.
	//Now we can load the stored values into the output variables.
	refpose_string=refpose_string_out;
	refpose_resnumber=refpose_resnumber_out;

	//Parse the offset:
	if ( hasoffset ) {
		std::istringstream ss( str.substr( endbracketposition+2, str.length() ) );
		ss >> refpose_offset;
		if ( str[endbracketposition+1] == '-' ) refpose_offset*=(-1);
	} else {
		refpose_offset=0;
	}

	//TR << "Parsing complete.  This string refers to ReferencePose " << refpose_string << ", residue " << refpose_resnumber << ", offset " << refpose_offset << ".  Returning true." << std::endl; //DELETE ME

	return true;
}

/// @brief Given the name of a ReferencePose object in the pose, a residue number in that reference pose, and a residue offset,
/// this function returns the Rosetta number of the corresponding residue in the pose.  Should throw an error if the ReferencePose
/// doesn't exist in the pose, or 0 if no corresponding residue exists in the pose.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
core::Size get_resnumber_from_reference_pose(
	std::string const &refpose_name,
	core::Size const refpose_number,
	signed long const refpose_offset,
	core::pose::Pose const &pose
) {
	signed long returnval( static_cast< signed long >( pose.corresponding_residue_in_current(refpose_number, refpose_name) ) + refpose_offset );
	if ( returnval <= 0 || returnval > static_cast< signed long >( pose.n_residue() ) ) return 0;
	return static_cast< core::Size >( returnval );
}

} //pose
} //core
