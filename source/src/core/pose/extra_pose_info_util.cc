// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/util.cc
/// @brief  Pose class utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Rhiju Das, Steven Lewis, Vikram K. Mulligan


// Unit header
#include <core/pose/extra_pose_info_util.hh>

// Package headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/PositionConservedResiduesStore.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/util.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/Exceptions.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/raw_data/DisulfideFile.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringIntegerMap.hh>
#include <basic/datacache/CacheableStringMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.string.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/string_constants.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/functional/hash.hpp>

namespace core {
namespace pose {

static basic::Tracer TR( "core.pose.extra_pose_info_util" );

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details Crude way to guess secondary structure given a pose. This function
/// sets the sec_struct array of pose.conformation_ to the result of the
/// guesswork. This has been ported directly from rosetta++.
void
set_ss_from_phipsi(
	pose::Pose & pose
) {
	// ss       :ss = 1 helix, ss = 2 sheet, ss = 3 other
	const int sstemp_offset=3;
	utility::vector1 < int > sstemp( sstemp_offset*2 + pose.size() );
	utility::vector1 < int > ss( pose.size() );

	sstemp[sstemp_offset-1] = 3; // assign loop to fictious residues at ends of chain
	sstemp[sstemp_offset+0] = 3;
	sstemp[sstemp_offset+pose.size()+1] = 3;
	sstemp[sstemp_offset+pose.size()+2] = 3;

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !pose.residue_type( i ).is_protein() ) { // make sure we don't inquire about the phi/psi of a non-protein residue
			sstemp[sstemp_offset+i] = 3;
		} else {
			if ( pose.phi(i) < -20.0 && pose.psi(i) > -90.0 && pose.psi(i) < -10.0 ) {
				sstemp[sstemp_offset+i] = 1;
			} else if ( pose.phi(i) < -20.0 && (pose.psi(i) > 20.0 || pose.psi(i) < -170.0) ) {
				sstemp[sstemp_offset+i] = 2;
			} else {
				sstemp[sstemp_offset+i] = 3;
			}
		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( sstemp[sstemp_offset+i] == 2 ) {
			if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i+1] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i-1] == 2 && sstemp[sstemp_offset+i-2] == 2 ) {
				ss[i] = 2;
			} else if ( sstemp[sstemp_offset+i+1] == 2 && sstemp[sstemp_offset+i+2] == 2 ) {
				ss[i] = 2;
			} else {
				ss[i] = 3;
			}
		} else if ( sstemp[sstemp_offset+i] == 1 ) {
			if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i+1] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i-1] == 1 && sstemp[sstemp_offset+i-2] == 1 ) {
				ss[i] = 1;
			} else if ( sstemp[sstemp_offset+i+1] == 1 && sstemp[sstemp_offset+i+2] == 1 ) {
				ss[i] = 1;
			} else {
				ss[i] = 3;
			}
		} else {
			ss[i] = 3;
		}
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( ss[i] == 1 ) {
			pose.set_secstruct(i,'H');
		} else if ( ss[i] == 2 ) {
			pose.set_secstruct(i,'E');
		} else {
			pose.set_secstruct(i,'L');
		}
	}
}

/// @brief  Reads the comments from the pdb file and adds it into comments
void
read_comment_pdb(
	std::string const &file_name,
	core::pose::Pose & pose
) {
	utility::io::izstream data(file_name);
	if ( !data ) {
		TR << "ERROR! PDB FILE '" << file_name << "' NOT FOUND!!!" <<std::endl;
		utility_exit_with_message("Cannot open PDB file.");
	}
	std::string line;
	while ( getline( data, line ) ) {
		if ( line != "##Begin comments##" ) continue;

		getline( data, line );
		while ( line != "##End comments##" ) {
			//TR<<"Testing read comments! :"<<line<<std::endl;
			utility::vector1<std::string> comment_line(utility::string_split(line,' '));
			core::pose::add_comment(pose,comment_line[1],comment_line[2]);
			getline( data, line );
		}
	}
}

void
dump_comment_pdb(
	std::string const &file_name,
	core::pose::Pose const& pose
) {
	std::ofstream out( file_name.c_str() );
	pose.dump_pdb(out);
	// verbose output
	out << "END\n";
	out << "##Begin comments##" << std::endl;
	using namespace std;
	map< string, string > const comments = core::pose::get_all_comments(pose);
	for ( auto const & comment : comments ) {
		out << comment.first<<" "<<comment.second << std::endl;
	}
	out << "##End comments##" << std::endl;
	out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
hasPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name )
{
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		return false;
	}

	CacheableStringFloatMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	debug_assert( data.get() != nullptr );

	return (  data->map().find( name ) != data->map().end() );
}

bool
hasPoseExtraScore_str(
	core::pose::Pose const & pose,
	std::string const & name )
{
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		return false;
	}

	CacheableStringMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
	debug_assert( data.get() != nullptr );

	return (  data->map().find( name ) != data->map().end() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	core::Real & value
) {
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		return false;
	}

	CacheableStringFloatMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
	debug_assert( data.get() != nullptr );

	auto it = data->map().find( name );
	if ( it == data->map().end() ) return false;

	value = it->second;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Real
getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name
) {
	Real value;
	runtime_assert( getPoseExtraScore( pose, name, value ) );
	return value;
}

bool getPoseExtraScore(
	core::pose::Pose const & pose,
	std::string const & name,
	std::string & value
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		return false;
	}

	CacheableStringMapCOP data
		= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
	debug_assert( data.get() != nullptr );

	auto it = data->map().find( name );
	if ( it == data->map().end() ) {
		return false;
	}
	value = it->second;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	core::Real value
) {
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringFloatMap() )

		);
	}

	CacheableStringFloatMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA) );

	runtime_assert( data.get() != nullptr );
	data->map()[name] = value;
}


void setPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name,
	std::string const & value
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has one of these.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )
		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA) );

	runtime_assert( data.get() != nullptr );
	data->map()[name] = value;
}

void add_comment(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );

	runtime_assert( data.get() != nullptr );
	data->map()[key] = val;
} // add_comment

void add_score_line_string(
	core::pose::Pose & pose,
	std::string const & key,
	std::string const & val
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	// make sure that the pose has a map of strings
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )
		);
	}

	CacheableStringMapOP data
		=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
		( pose.data().get_ptr(core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS) );

	runtime_assert( data.get() != nullptr );
	data->map()[key] = val;
}

void clearPoseExtraScores(
	core::pose::Pose & pose
) {

	using basic::datacache::DataCache_CacheableData;

	{
		using basic::datacache::CacheableStringFloatMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringFloatMap() )
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )
		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::STRING_MAP,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}

	{
		using basic::datacache::CacheableStringMap;
		pose.data().set(
			core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS,
			DataCache_CacheableData::DataOP( new CacheableStringMap() )

		);
	}
}

void clearPoseExtraScore(
	core::pose::Pose & pose,
	std::string const & name
) {
	using basic::datacache::CacheableStringFloatMap;
	using basic::datacache::CacheableStringFloatMapOP;
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;
	using basic::datacache::DataCache_CacheableData;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		CacheableStringFloatMapOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringFloatMap >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		debug_assert( data.get() != nullptr );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}

	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		CacheableStringMapOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		debug_assert( data.get() != nullptr );

		if ( data->map().find( name ) != data->map().end() ) data->map().erase( name );
	}
}

bool get_comment(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > comment_map = get_all_comments( pose );

	if ( comment_map.find( key ) == comment_map.end() ) {
		return false;
	}

	val = comment_map[ key ];
	return true;
}

bool get_score_line_string(
	core::pose::Pose const & pose,
	std::string const & key,
	std::string & val
) {
	using std::map;
	using std::string;
	map< string, string > score_line_strings_map = get_all_score_line_strings( pose );

	if ( score_line_strings_map.find( key ) == score_line_strings_map.end() ) {
		return false;
	}

	val = score_line_strings_map[ key ];
	return true;
}

void delete_comment(
	core::pose::Pose & pose,
	std::string const & key
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapOP;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapOP data
			=  utility::pointer::dynamic_pointer_cast< CacheableStringMap >
			( pose.data().get_ptr(core::pose::datacache::CacheableDataType::STRING_MAP) );
		std::map< std::string, std::string >::iterator it;
		it = data->map().find(key);
		if ( it != data->map().end() ) {
			data->map().erase(it);
		}
	}
}

std::map< std::string, std::string >
get_all_score_line_strings(
	core::pose::Pose const & pose
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > score_line_strings;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) ) {
		CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::SCORE_LINE_STRINGS ) );
		score_line_strings = data->map();
		runtime_assert( data.get() != nullptr );
	}
	return score_line_strings;
}

std::map< std::string, std::string >
get_all_comments(
	core::pose::Pose const & pose
) {
	using basic::datacache::CacheableStringMap;
	using basic::datacache::CacheableStringMapCOP;

	std::map< std::string, std::string > comments;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::STRING_MAP ) ) {
		CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STRING_MAP ) );
		comments = data->map();
		runtime_assert( data.get() != nullptr );
	}
	return comments;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose ) {
	using namespace basic::options;

	std::string filename( option[ OptionKeys::in::file::psipred_ss2 ]().name() );
	return read_psipred_ss2_file( pose, filename );
}

utility::vector1< char > read_psipred_ss2_file( pose::Pose const & pose, std::string const & filename ) {
	using namespace basic::options;

	utility::vector1< char > secstructs;
	utility::io::izstream data( filename );

	if ( !data ) {
		TR.Warning << "Cannot open psipred_ss2 file " << filename << std::endl;
		return secstructs;
	}

	std::string line;
	Size count(0);
	while ( getline( data, line ) ) {
		if ( line[0] == '#' || line == "" ) {
			continue;
		}
		std::istringstream line_stream( line );
		Size pos;
		char aa, sec;
		line_stream >> pos >> aa >> sec;
		count++;
		if ( sec != 'H' && sec != 'E' && sec != 'C' ) {
			TR.Warning << "unrecognized secstruct char : " << sec << " at seqpos " << count << std::endl;
		}
		if ( sec == 'C' ) {
			secstructs.push_back( 'L' );
		} else {
			secstructs.push_back( sec );
		}
	}

	// chu get number of protein residues
	Size nres=0;
	for ( Size i =1 ; i <= pose.size(); ++i ) {
		if ( pose.residue(i).is_protein() ) nres++;
	}

	debug_assert( secstructs.size() == nres);
	if ( secstructs.size() != nres ) {
		secstructs.clear();
	}

	return secstructs;
}

std::string tag_from_pose( core::pose::Pose const & pose ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;

	std::string tag( "empty_tag" );
	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag =
			static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}

	return tag;
}

std::string extract_tag_from_pose( core::pose::Pose &pose )
{
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	using basic::datacache::CacheableStringOP;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		CacheableStringOP data =  utility::pointer::dynamic_pointer_cast< CacheableString > (  (pose.data().get_ptr( ( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG  )) ));
		if ( data.get() == nullptr ) return std::string("UnknownTag");
		else                      return data->str();
	}

	return std::string("UnknownTag");
}


void tag_into_pose( core::pose::Pose & pose, std::string const & tag ) {
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;
	using basic::datacache::CacheableString;
	using basic::datacache::DataCache_CacheableData;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new CacheableString(tag) )
	);
}

///////////////////////////////////////////////////////////////////////////////
void
set_output_res_and_chain( pose::Pose & extended_pose,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & output_resnum_and_chain )
{
	using namespace pose;
	utility::vector1< int > const & output_res_num = std::get< 0 >( output_resnum_and_chain );
	utility::vector1< char > const & output_chain = std::get< 1 >( output_resnum_and_chain );
	utility::vector1< std::string > const & output_segid = std::get< 2 >( output_resnum_and_chain );
	if ( output_res_num.size() == 0 ) return;
	runtime_assert( output_res_num.size() == extended_pose.size() );

	PDBInfoOP pdb_info( new PDBInfo( extended_pose ) );
	pdb_info->set_numbering( output_res_num );

	bool chain_interesting( false );
	for ( Size n = 1; n <= output_chain.size(); n++ ) {
		if ( output_chain[ n ] != ' ' ) chain_interesting = true;
	}
	if ( chain_interesting ) pdb_info->set_chains( output_chain );

	bool segid_interesting( false );
	for ( Size n = 1; n <= output_segid.size(); n++ ) {
		if ( output_segid[ n ] != "    " ) segid_interesting = true;
	}
	if ( segid_interesting ) pdb_info->set_segment_ids( output_segid );

	extended_pose.pdb_info( pdb_info );
}


} // pose
} // core
