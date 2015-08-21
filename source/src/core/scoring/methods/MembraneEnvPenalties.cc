// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RMS_Energy.cc
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/MembraneEnvPenalties.hh>
#include <core/scoring/methods/MembraneEnvPenaltiesCreator.hh>

// Package headers

#include <core/scoring/MembranePotential.hh>
#include <core/scoring/ScoringManager.hh>
//#include <core/scoring/rms_util.hh>

//#include <core/io/pdb/pose_io.hh>

//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>


#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


//#include <basic/prof.hh>
//#include <utility/exit.hh>

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MembraneEnvPenalties class,
/// never an instance already in use
methods::EnergyMethodOP
MembraneEnvPenaltiesCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MembraneEnvPenalties );
}

ScoreTypes
MembraneEnvPenaltiesCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( Menv_non_helix );
	sts.push_back( Menv_termini );
	sts.push_back( Menv_tm_proj );
	return sts;
}


/// c-tor
MembraneEnvPenalties::MembraneEnvPenalties() :
	parent( EnergyMethodCreatorOP( new MembraneEnvPenaltiesCreator ) ),
	potential_( ScoringManager::get_instance()->get_MembranePotential() )
{}


/// clone
EnergyMethodOP
MembraneEnvPenalties::clone() const
{
	return EnergyMethodOP( new MembraneEnvPenalties() );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the RMS difference between native_pose_ (provided by
/// the option -in::file::native and the given Pose. The actual energy calculation
/// is the difference between the RMSD and the target RMSD. Target RMSD is specified
/// the option -score::rms_target.


/*
void
MembraneEnvPenalties::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
// compute interpolated number of neighbors at various distance cutoffs
pose.update_residue_neighbors();
potential_.compute_centroid_environment( pose );
potential_.compute_membrane_embedding( pose );

}
*/
void
MembraneEnvPenalties::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	if ( potential_.Menv_penalties() ) { //bw quick hack before putting them as individual scoring terms....
		Real tm_projection(0);
		Real non_helix_pen(0);
		Real termini_pen(0);
		potential_.tm_projection_penalty(pose,tm_projection);
		potential_.non_helix_in_membrane_penalty(pose, non_helix_pen);
		potential_.termini_penalty(pose,termini_pen);
		emap[ Menv_non_helix ]=non_helix_pen;
		emap[ Menv_termini ]=termini_pen;
		emap[ Menv_tm_proj ]=tm_projection;


		//  std::cout << "Menv_penalties (tm_projection+hbond_pen+termini_pen+10) " << tm_projection << " " << hbond_pen << " " << termini_pen << std::endl;
	}
	potential_.finalize( pose );
	// totals[ rms ]  = std::abs( rms_target_ - rmsd );

	// PROF_STOP( basic::RMS );
}
core::Size
MembraneEnvPenalties::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
