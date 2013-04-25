// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RG_Energy.cc
/// @brief  Radius of gyration energy function definition. Returns -1 * RG for a given Pose.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/RG_Energy.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <core/scoring/EnvPairPotential.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


// Utility headers
#include <basic/prof.hh>
#include <ObjexxFCL/format.hh>


// C++


namespace core {
namespace scoring {
namespace methods {

/// THIS CLASS IS DEPRICATED.
/// USE RG_Energy_Fast
RG_Energy::RG_Energy()
{
	//add_score_type( rg );
}


/// clone
EnergyMethodOP
RG_Energy::clone() const
{
	return new RG_Energy();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the radius of gyration and place the answer into
/// totals[ rg ].
void
RG_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	PROF_START( basic::RG );

	Size const nres( pose.total_residue() );

	///////////////////////////////////////
	//
	// RG SCORE
	core::Real rg_score = 0;
	for( Size i = 1; i <= nres; ++i ) { // outer loop over residues
		Vector const v( pose.residue(i).nbr_atom_xyz() );
		for ( Size j = i + 1; j <= nres; ++j ) {
			rg_score += distance_squared( v, pose.residue(j).nbr_atom_xyz() );
		}
	}

	rg_score /= ( pose.total_residue() * pose.total_residue() - pose.total_residue() );
	totals[ rg ] = sqrt(rg_score);

	PROF_STOP( basic::RG );
}

}
}
}
