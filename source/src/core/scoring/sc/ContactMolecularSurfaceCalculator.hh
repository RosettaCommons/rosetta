// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/sc/ContactMolecularSurfaceCalculator.hh
/// @brief  Headers for the Weighted Contact Molecular Surface Calculator
/// @author Longxing Cao (longxing@uw.edu)

#ifndef INCLUDED_core_scoring_sc_ContactMolecularSurfaceCalculator_hh
#define INCLUDED_core_scoring_sc_ContactMolecularSurfaceCalculator_hh

//
// Note: This almost use the same code as ShapeComplementarityCalculator
// and it returns the weighted contact molecular surface of the target.
// Just avoid of messing up the ShapeComplementarityCalculator, so just
// copied most of the code.


// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/sc/ContactMolecularSurfaceCalculator.fwd.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>

//// C++ headers
#include <vector>
#include <string>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace sc {


////////////////////////////////////////////////////////////
// Contact Molecular Surface  Calculator class definition
////////////////////////////////////////////////////////////

class ContactMolecularSurfaceCalculator : public MolecularSurfaceCalculator {

public:

	ContactMolecularSurfaceCalculator();
	~ContactMolecularSurfaceCalculator() override;

	MolecularSurfaceCalculator::ScValue CalcContactArea();

protected:
	// Surface generation configuration
	ScValue CalcContactMolecularSurface( std::vector<DOT> const &, std::vector<DOT> const & );
	int AssignAttentionNumbers( std::vector<Atom>& atom ) override;

private:

	DOT const *CalcNeighborDistanceFindClosestNeighbor( DOT const & dot1, std::vector<const DOT*> const & their_dots );


};

} //namespace sc
} //namespace filters
} //namespace protocols

#endif

