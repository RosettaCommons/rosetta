// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ContactMolecularSurfaceCalculator.cc
/// @brief    Headers for the Shape Complementarity Calculator
/// @details  Lawrence & Coleman shape complementarity calculator (based on CCP4's sc)
/// @author   Longxing Cao <longxing@uw.edu>

#ifndef INCLUDED_core_scoring_sc_ContactMolecularSurfaceCalculator_cc
#define INCLUDED_core_scoring_sc_ContactMolecularSurfaceCalculator_cc

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ContactMolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>


// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <utility/excn/Exceptions.hh>

static basic::Tracer TR( "core.scoring.sc.ContactMolecularSurfaceCalculator" );

using namespace core;

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////////////////////
// Public class functions
////////////////////////////////////////////////////////////////////////////

/// @brief
/// ContactMolecularSurfaceCalculator constructor, initializes default settings
ContactMolecularSurfaceCalculator::ContactMolecularSurfaceCalculator() :
	MolecularSurfaceCalculator()
{
}

ContactMolecularSurfaceCalculator::~ContactMolecularSurfaceCalculator() = default;


MolecularSurfaceCalculator::ScValue ContactMolecularSurfaceCalculator::CalcContactArea()
{
	try {
		run_.results.valid = 0;

		if ( run_.atoms.empty() ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined");
		}
		if ( !run_.results.surface[0].nAtoms ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined for molecule 1");
		}
		if ( !run_.results.surface[1].nAtoms ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined for molecule 2");
		}

		// Determine and assign the attention numbers for each atom
		AssignAttentionNumbers(run_.atoms);

		GenerateMolecularSurfaces();

		if ( !run_.dots[0].size() || !run_.dots[1].size() ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No molecular dots generated!");
		}

		ScValue area = CalcContactMolecularSurface( run_.dots[0], run_.dots[1] );

		return area;
	} catch ( ShapeComplementarityCalculatorException & e ) {
		TR.Error << "Failed: " << e.error << std::endl;
	}

	return -1;
}

////////////////////////////////////////////////////////////////////////////
// Protected class functions
////////////////////////////////////////////////////////////////////////////

// Determine assign the attention numbers for each atom
int ContactMolecularSurfaceCalculator::AssignAttentionNumbers(std::vector<Atom> & )
{
	std::vector<Atom>::iterator pAtom1, pAtom2;

	for ( pAtom1 = run_.atoms.begin(); pAtom1 < run_.atoms.end(); ++pAtom1 ) {
		// find nearest neighbour in other molecule
		ScValue dist_min = 99999.0, r;
		for ( pAtom2 = run_.atoms.begin(); pAtom2 < run_.atoms.end(); ++pAtom2 ) {
			if ( pAtom1->molecule == pAtom2->molecule ) {
				continue;
			}
			r = pAtom1->distance(*pAtom2);
			if ( r < dist_min ) {
				dist_min = r;
			}
		}

		// check if within separator distance
		if ( dist_min >= settings.sep ) {
			// TR.Debug << "Atom ATTEN_BLOCKER: " << pAtom1->natom << std::endl;
			// too _far_ away from other molecule, blocker atom only
			pAtom1->atten = ATTEN_BLOCKER;
			++run_.results.surface[pAtom1->molecule].nBlockedAtoms;
		} else {
			// potential interface or neighbouring atom
			pAtom1->atten = ATTEN_BURIED_FLAGGED;
			++run_.results.surface[pAtom1->molecule].nBuriedAtoms;
		}
	}

	return 1;
}


ContactMolecularSurfaceCalculator::ScValue ContactMolecularSurfaceCalculator::CalcContactMolecularSurface(
	std::vector<DOT> const & my_dots,
	std::vector<DOT> const & their_dots)
{
	ScValue area = 0.0;

	if ( my_dots.empty() ) return 0.0;
	std::vector<DOT const *> buried_their_dots;
	for ( auto idot = their_dots.begin(); idot < their_dots.end(); ++idot ) {
		DOT const &dot = *idot;
		if ( dot.buried ) buried_their_dots.push_back( &dot );
	}

	for ( auto idot = my_dots.begin(); idot < my_dots.end(); ++idot ) {
		DOT const &dot = *idot;
		if ( !dot.buried ) continue;
		DOT const *neighbor = nullptr;
		neighbor = CalcNeighborDistanceFindClosestNeighbor(dot, buried_their_dots);
		if ( !neighbor ) continue;
		ScValue distmin = neighbor->coor.distance(dot.coor);
		area += dot.area * exp( - pow(distmin, 2) * settings.weight );
	}
	return area;
}

// Find closest neighbor dot for a given dot
DOT const * ContactMolecularSurfaceCalculator::CalcNeighborDistanceFindClosestNeighbor(
	DOT const & dot1,
	std::vector<DOT const*> const & their_dots
) {
	MolecularSurfaceCalculator::ScValue distmin = 999999.0, d;
	DOT const *neighbor = nullptr;

	// Loop over the entire surface: find and flag neighbour of each point
	// that we're interested in and store nearest neighbour pointer

	for ( auto idot2 = their_dots.begin();
			idot2 < their_dots.end(); ++idot2 ) {
		DOT const &dot2 = **idot2;
		if ( !dot2.buried ) continue;
		d = dot2.coor.distance_squared(dot1.coor);
		if ( d <= distmin ) {
			distmin = d;
			neighbor = *idot2;
		}
	}
	return neighbor;
}


} // namespace sc
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_sc_ContactMolecularSurfaceCalculator_cc

// END //
