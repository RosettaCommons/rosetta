// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SurfaceEnergy.cc
/// @author Ron Jacak


// Unit headers
#include <core/pack/interaction_graph/SurfaceEnergy.hh>
#include <core/pack/interaction_graph/SurfaceEnergyCreator.hh>

// Package headers
#include <core/pack/interaction_graph/SurfacePotential.hh>
//#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "core.pack.interaction_graph.SurfaceEnergy" );

// C++ headers


namespace core {
namespace pack {
namespace interaction_graph {


/// @details This must return a fresh instance of the SurfaceEnergy class,
/// never an instance already in use
scoring::methods::EnergyMethodOP
SurfaceEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const &
) const {
	return new SurfaceEnergy;
}

scoring::ScoreTypes
SurfaceEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( scoring::surface );
	return sts;
}


SurfaceEnergy::SurfaceEnergy() :
	parent( new SurfaceEnergyCreator )
{}


scoring::methods::EnergyMethodOP
SurfaceEnergy::clone() const {
	return new SurfaceEnergy();
}


void
SurfaceEnergy::setup_for_scoring( pose::Pose & /*pose*/, scoring::ScoreFunction const & /*sf*/ ) const {
	// since this is a fake EnergyMethod, don't do anything here
	// is it even necessary to implement this method as empty?
}


void
SurfaceEnergy::residue_energy( conformation::Residue const & /*rsd*/, pose::Pose const & /*pose*/, scoring::EnergyMap & /*emap*/ ) const {

	// if this were a real term, the code here might look like the following
	//Real surface_score( 0.0 );
	//surface_score = evaluate_env_score( pose, rsd )
	//emap [ surface ] += surface_score

}

// SurfaceEnergy is non-pairwise decomposable, so it can only be calculated all at once, not on a residue by residue
// basis. For packing, it uses a specialized InteractionGraph that updates the score efficiently for substitutions.
void
SurfaceEnergy::finalize_total_energy( pose::Pose & pose, scoring::ScoreFunction const &, scoring::EnergyMap & totals ) const {

	// don't run if minimizing, non-differentiable
	if ( ! pose.energies().use_nblist() ) {
		core::Real surface_score;
		SurfacePotential::get_instance()->compute_pose_surface_energy( pose, surface_score );
		//TR << "surface score: " << surface_score << std::endl;
		totals[ scoring::surface ] = surface_score;
	}

}
core::Size
SurfaceEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace interaction_graph
} // pack
} // core
