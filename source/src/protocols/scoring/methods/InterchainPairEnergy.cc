// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/InterchainPairEnergy.cc
/// @brief  Statistically derived rotamer pair potentials
/// @detailed For docking (or between chains) only those residues at the interface
///						and between the two interfaces need to be evaluated
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/methods/InterchainPairEnergy.hh>
#include <protocols/scoring/methods/InterchainPairEnergyCreator.hh>

// Package headers
#include <protocols/scoring/InterchainPotential.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



// Utility headers



// C++

namespace protocols {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the InterchainPairEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
InterchainPairEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return new InterchainPairEnergy;
}

core::scoring::ScoreTypes
InterchainPairEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( interchain_pair );
	sts.push_back( interchain_vdw );
	return sts;
}


/// c-tor
InterchainPairEnergy::InterchainPairEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new InterchainPairEnergyCreator ) ),
	potential_( InterchainPotential::get_instance() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
InterchainPairEnergy::clone() const
{
	return new InterchainPairEnergy();
}

///
void
InterchainPairEnergy::setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	potential_->compute_centroid_environment( pose );
	potential_->compute_interface( pose );
	pose.update_residue_neighbors();
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
/// @brief calculate the pair_score between chains
void
InterchainPairEnergy::residue_pair_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core;
	using namespace core::scoring;
	Real pair_score ( 0.0 ), vdw_score ( 0.0 );

	/// this is probably really slow, interface should probably be kept track of in the interaction graph
	potential_->evaluate_pair_and_vdw_score( pose, rsd1, rsd2, pair_score, vdw_score );

	/// divide by two to account for double counting that will occur
	emap[ interchain_pair ] += pair_score;
	emap[ interchain_vdw ] += vdw_score;
}

// is there a better way to get the contact score??
void
InterchainPairEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &// emap
) const
{
	// sets calculated from the CenPairInfo to false
	potential_->finalize( pose );
}

/// @brief Energy distance cutoff
core::Distance
InterchainPairEnergy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 6.0 from cutoffs in centroid params files
// 	return 0.0; /// since all the cutoffs for centroid mode are rolled into the cendist check
}
core::Size
InterchainPairEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
