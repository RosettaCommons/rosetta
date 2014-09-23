// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HPatchEnergy.cc
/// @author Ron Jacak

// Unit headers
#include <core/pack/interaction_graph/HPatchEnergy.hh>
#include <core/pack/interaction_graph/HPatchEnergyCreator.hh>

// Package headers
#include <core/pack/interaction_graph/SurfacePotential.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "core.pack.interaction_graph.HPatchEnergy" );

// C++ headers


namespace core {
namespace pack {
namespace interaction_graph {

/// @details This must return a fresh instance of the HPatchEnergy class, never an instance already in use
scoring::methods::EnergyMethodOP
HPatchEnergyCreator::create_energy_method( scoring::methods::EnergyMethodOptions const & ) const {
	return new HPatchEnergy;
}

scoring::ScoreTypes
HPatchEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( scoring::hpatch );
	return sts;
}


HPatchEnergy::HPatchEnergy() :
	parent( scoring::methods::EnergyMethodCreatorOP( new HPatchEnergyCreator ) )
{}


scoring::methods::EnergyMethodOP
HPatchEnergy::clone() const {
	return new HPatchEnergy();
}


void
HPatchEnergy::setup_for_scoring( pose::Pose & /*pose*/, scoring::ScoreFunction const & /*sf*/ ) const {
	// since this is a fake EnergyMethod, don't do anything here
	// is it even necessary to implement this method as empty?
}


void
HPatchEnergy::residue_energy( conformation::Residue const & /*rsd*/, pose::Pose const & /*pose*/, scoring::EnergyMap & /*emap*/ ) const {

	// if this were a real term, the code here might look like the following
	//Real hpatch_score( 0.0 );
	//hpatch_score = evaluate_score( pose, rsd )
	//emap [ hpatch ] += hpatch_score

}

// HPatchEnergy is non-pairwise decomposable, so it can only be calculated all at once, not on a residue by residue
// basis. For packing, it uses a specialized InteractionGraph that updates the score efficiently for substitutions.
void
HPatchEnergy::finalize_total_energy( pose::Pose & pose, scoring::ScoreFunction const &, scoring::EnergyMap & totals ) const {

	// don't run if minimizing, non-differentiable
	if ( ! pose.energies().use_nblist() ) {
		core::Real hpatch_score = 0.0;
		std::map< Size, std::pair< Real, Real > > patch_scores;
		std::map< Size, utility::vector1< id::AtomID > > atoms_in_patches;
		SurfacePotential::get_instance()->compute_pose_hpatch_score( pose, hpatch_score, patch_scores, atoms_in_patches );
		//TR << "hpatch score: " << hpatch_score << std::endl;
		totals[ scoring::hpatch ] = hpatch_score;
	}

}
core::Size
HPatchEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace interaction_graph
} // namespace pack
} // namespace core
