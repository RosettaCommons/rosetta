// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/io/util.cc
/// @brief Util functions for Input and Output.  Very general IO should go to utility/io.
///   These should be related to core in a deep way.
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

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>

namespace core {
namespace io {

using namespace ObjexxFCL::format; // AUTO USING NS

core::Real restrict_prec( core::Real inval )
{
	if ( inval >= 1 || inval <= -1 ) { // Don't alter value, as the default precision of 6 works fine, and we avoid rounding artifacts
		return inval;
	}
	core::Real outval;
	std::stringstream temp;
	temp << std::fixed << std::setprecision(5) << inval;
	temp >> outval;
	return outval;
}

/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string extract_scores(
	core::pose::Pose const & pose,
	std::string const &filename
) {
	std::stringstream out;
	extract_scores(pose, out, filename);
	return out.str();
}

void extract_scores(
	core::pose::Pose const & pose,
	std::stringstream & out,
	std::string const &filename
)
{
	if ( basic::options::option[ basic::options::OptionKeys::out::file::no_scores_in_pdb ] ) {
		return;
	}

	//This is shamelessly refactored from the older JobDistributor; Jobdistributors.hh:1018; SVN 25940
	// APL: Moving this job-independent code into a central location.
	// Which score terms to use
	core::scoring::EnergyMap weights = pose.energies().weights();
	typedef utility::vector1<core::scoring::ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if ( weights[ii] != 0 ) score_types.push_back(ii);
	}
	// This version is formatted for easy parsing by R, Excel, etc.
	out << "# All scores below are weighted scores, not raw scores.\n";
	out << "#BEGIN_POSE_ENERGIES_TABLE " << filename << std::endl;
	out << "label";
	BOOST_FOREACH ( core::scoring::ScoreType score_type, score_types ) {
		out << " " << name_from_score_type(score_type);
	}
	out << " total\n";
	out << "weights";
	BOOST_FOREACH ( core::scoring::ScoreType score_type, score_types ) {
		out << " " << weights[score_type];
	}
	out << " NA\n";
	out << "pose";
	core::Real pose_total = 0.0;
	if ( pose.energies().energies_updated() ) {
		BOOST_FOREACH ( core::scoring::ScoreType score_type, score_types ) {
			core::Real score = (weights[score_type] * pose.energies().total_energies()[ score_type ]);
			out << " " << restrict_prec(score);
			pose_total += score;
		}
		out << " " << restrict_prec(pose_total) << "\n";
		for ( core::Size j = 1, end_j = pose.total_residue(); j <= end_j; ++j ) {
			core::Real rsd_total = 0.0;
			out << pose.residue(j).name() << "_" << j;
			BOOST_FOREACH ( core::scoring::ScoreType score_type, score_types ) {
				core::Real score = (weights[score_type] * pose.energies().residue_total_energies(j)[ score_type ]);
				out << " " << restrict_prec(score);
				rsd_total += score;
			}
			out << " " << restrict_prec(rsd_total) << "\n";
		}
	}
	out << "#END_POSE_ENERGIES_TABLE " << filename << std::endl;
}


/// @brief Write Pose energies information into a string and return it.
/// @details Added during the 2016 Chemical XRW.  This is a bloody mess.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string extract_extra_scores(
	pose::Pose const & pose
) {
	std::stringstream out;
	extract_extra_scores(pose, out);
	return out.str();
}

void extract_extra_scores(
	pose::Pose const & pose,
	std::stringstream & out
)
{
	// ARBITRARY_STRING_DATA
	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) ) {
		basic::datacache::CacheableStringMapCOP data
			= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_STRING_DATA ) );
		assert( data.get() != NULL );

		for ( std::map< std::string, std::string >::const_iterator it( data->map().begin() ), end( data->map().end() );
				it != end;
				++it ) {
			//TR << it->first << " " << it->second << std::endl;
			out << it->first << " " << it->second << std::endl;
		}
	}

	// ARBITRARY_FLOAT_DATA
	if ( pose.data().has( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) {
		basic::datacache::CacheableStringFloatMapCOP data
			= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringFloatMap const >
			( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) );
		assert( data.get() != NULL );

		for ( std::map< std::string, float >::const_iterator it( data->map().begin() ), end( data->map().end() );
				it != end;
				++it ) {
			//TR << it->first << " " << it->second << std::endl;
			out << it->first << " " << it->second << std::endl;
		}
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
