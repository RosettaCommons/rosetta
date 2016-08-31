// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/util.cc
/// @brief Util functions for Input and Output.  Very general IO should go to utility/io.
///   These should be related to core in a deep way or not able to be called from utility.
///
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

// Unit headers
#include <core/io/util.hh>

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>


// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/datacache/CacheableStringFloatMap.hh>
#include <basic/datacache/CacheableStringMap.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>

namespace core {
namespace io {

using namespace ObjexxFCL::format; // AUTO USING NS



/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu) Jared Adolf-Bryfogle (jadolfbr@gmail.com)
std::string pose_energies_from_sfr(
	StructFileRep const & sfr
) {
	std::stringstream out;
	pose_energies_from_sfr(sfr, out);
	return out.str();
}

void pose_energies_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
)
{
	using namespace core::io::pose_to_sfr;

	// This version is formatted for easy parsing by R, Excel, etc.

	utility::vector1< std::string > const & score_names = sfr.score_table_labels();
	utility::vector1< std::vector< std::string > > const & score_lines = sfr.score_table_lines();
	
	if (score_names.size() == 0 || score_lines.size() == 0) return; //Was not extracted!
	
	out << "# All scores below are weighted scores, not raw scores.\n";
	
	if (! sfr.score_table_filename().empty()){
		out << "#BEGIN_POSE_ENERGIES_TABLE " << sfr.score_table_filename() << std::endl;
	}
	else {
		out << "#BEGIN_POSE_ENERGIES_TABLE " << std::endl;
	}
	
	out << "label";
	
	BOOST_FOREACH ( std::string score_name, score_names ) {
		out << " " << score_name;
	}
	out << "\n";
	out << "weights";
	utility::vector1< core::Real > const & score_weights = sfr.score_table_weights();
	
	BOOST_FOREACH ( core::Real weight, score_weights ) {
		out << " " << weight;
	}
	out << " NA\n";
	
	
	BOOST_FOREACH( std::vector<std::string> score_line, score_lines){
		std::string line = "";
		BOOST_FOREACH( std::string column, score_line){
			line = line+" "+column;
		}
		line = utility::strip(line);
		out << line << "\n";
	}
	if (! sfr.score_table_filename().empty()){
		out << "#END_POSE_ENERGIES_TABLE " << sfr.score_table_filename() << std::endl;
	}
	else {
		out << "#END_POSE_ENERGIES_TABLE " << std::endl;
	}
}


/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu) + Jared Adolf-Bryfogle (jadolfbr@gmail.com)
std::string pose_data_cache_from_sfr(
	StructFileRep const & sfr
) {
	std::stringstream out;
	pose_data_cache_from_sfr(sfr, out);
	return out.str();
}

void pose_data_cache_from_sfr(
	StructFileRep const & sfr,
	std::stringstream & out
)
{

	//If either of these are empty, will not do anything.
	std::map< std::string, std::string > const & string_data = sfr.pose_cache_string_data();
	std::map< std::string,     float   > const & float_data =  sfr.pose_cache_float_data();
	
	// ARBITRARY_STRING_DATA
	for ( std::map< std::string, std::string >::const_iterator it( string_data.begin() ), end( string_data.end() );
				it != end;
				++it ){
		//TR << it->first << " " << it->second << std::endl;
		out << it->first << " " << it->second << std::endl;
	}

	// ARBITRARY_FLOAT_DATA
	for ( std::map< std::string, float >::const_iterator it( float_data.begin() ), end( float_data.end() );
				it != end;
				++it ) {
		//TR << it->first << " " << it->second << std::endl;
		out << it->first << " " << it->second << std::endl;
	}
}

void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices
){
	StructFileRepOptions options;
	pose_from_pose( new_pose, old_pose, residue_indices, options );
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
){
	using namespace chemical;
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	pose_from_pose( new_pose, old_pose, *residue_set,  residue_indices, options);
}


void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices
){
	StructFileRepOptions options;
	pose_from_pose( new_pose, old_pose, residue_set, residue_indices, options );
}


/// Creates a subpose from a pose, to include only certain
/// residues, using StructFileRep::init_from_pose() to construct the
/// pose, and build_pose_as_is1() to construct the pose
/// with the given options.
void
pose_from_pose(
	pose::Pose & new_pose,
	pose::Pose const & old_pose,
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< core::Size > const & residue_indices,
	StructFileRepOptions const & options
){
	core::io::pose_to_sfr::PoseToStructFileRepConverter converter;
	converter.init_from_pose( old_pose, residue_indices );
	pose_from_sfr::PoseFromSFRBuilder builder( residue_set.get_self_ptr(), options );
	builder.build_pose( *converter.sfr(), new_pose );

}



} //core
} //io
