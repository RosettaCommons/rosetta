// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/SurfacePotential.hh
/// @brief  Class which keeps reads the residue hydrophobic ASA database file and calculates surface residue energies.
/// @author Ron Jacak

#ifndef INCLUDED_core_pack_interaction_graph_SurfacePotential_hh
#define INCLUDED_core_pack_interaction_graph_SurfacePotential_hh

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

namespace core {
namespace pack {
namespace interaction_graph {

/// @brief With the traditional scoring hierarchy, classes like this one are created and accessed via the ScoringManager, which
/// is itself a Singleton class.  These "potential" classes are only created and initialized when the use of the EnergyMethod
/// these classes correspond is encountered.  No point in reading database files for a term if that term is not being used
/// in some score function.  However, the surface energy is used when users specify they want to use it on the command
/// line - NOT via a score function.  The score/energy is done within an interaction graph. One might ask why I just don't
/// put the logic for reading in the database file to the interaction graph init methods.  However, there will be cases
/// where I will want to just score a protein (and not do any design) where I will want the database file to be read in.
/// Scoring doesn't use interaction graphs, so if the code for that was located there, these values would not be read in.
/// Instead, I've decided to implement this as its own separate class.  It uses the Singleton design pattern so the database
/// will only get read in once during a run.
class SurfacePotential : public utility::SingletonBase< SurfacePotential >
{
public:
	friend class utility::SingletonBase< SurfacePotential >;

public:
	Real average_residue_hASA( chemical::AA aa_type, Size num_nbs );
	Real hASA_patch_energy( Real patch_area, Size num_nbs );

	void compute_residue_surface_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		scoring::EnergyMap & emap,
		Size resid,
		utility::vector1< Size > num_neighbors_
	);

	void compute_pose_surface_energy(
		pose::Pose const & pose,
		Real & surface_energy_
	);

	void compute_pose_surface_energy(
		pose::Pose const & pose,
		Real & total_surface_energy_,
		utility::vector1< Real > & residue_surface_energy_
	);

	core::Real hpatch_score( core::Real patch_area );

	/// @brief return the hpatch score for an entire pose
	Real compute_pose_hpatch_score(
		pose::Pose const & pose
	);

	void compute_pose_hpatch_score(
		pose::Pose const & pose,
		Real & total_hpatch_energy_,
		std::map< core::Size, std::pair< core::Real, core::Real > > & patch_scores_,
		std::map< core::Size, utility::vector1< id::AtomID > > & atoms_in_patches_
	);

	static const core::Size MAX_HPATCH_AREA;
	static const core::Real MAX_HPATCH_SCORE;
	static const core::Size HPATCH_SCORE_BIN_SIZE;

private:
	/// @brief private constructor
	SurfacePotential();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SurfacePotential * create_singleton_instance();

	void read_average_res_hASA_database_file();
	void read_hASA_score_database_file();
	void read_hpatch_score_database_file();

private:

	// outer vector holds AA's; inner vector holds neighbor counts. Residues always have at least 1 neighbor because
	// num_neighbors_counting_self() is always used to determine number of neighbors
	utility::vector1< utility::vector1 < Real > > res_to_average_hASA_;
	std::vector< utility::vector1 < Real > > hASA_to_score_;

	// vector which holds the values read from scoring/score_functions/SurfacePotential/hpatch_score.txt
	std::vector< core::Real > patcharea_to_score_;

	static const core::Size MAX_PATCH_SURFACE_AREA;
	static const core::Real MAX_SURFACE_ENERGY;
	static const core::Size SURFACE_EXPOSED_CUTOFF;
	static const core::Real INTERACTION_RADIUS;
	static const core::Size SURFACE_SCORE_BIN_SIZE;
	static const core::Size BURIED_RESIDUE_NO_HSASA_CUTOFF;

};

} // namespace interaction_graph
} // namespace pack
} // namespace core


#endif
