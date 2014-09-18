// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ./src/protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.cc
/// @brief  calss for helix-pairing potential
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit header
#include <protocols/fldsgn/potentials/sspot/NatbiasHelixPairPotential.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <numeric/numeric.functions.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.fldsgn.potentials.sspot.NatbiasHelixPairPotential", basic::t_info );

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

/// @brief default constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential():
	twist_( 15.0 )
{}


/// @brief value constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential( HelixPairingSetOP const hpairset ):
	twist_( 15.0 )
{}


/// @brief copy constructor
NatbiasHelixPairPotential::NatbiasHelixPairPotential( NatbiasHelixPairPotential const & src ):
	ReferenceCount(),
	twist_( src.twist_ )
{}


/// @brief default destructor
NatbiasHelixPairPotential::~NatbiasHelixPairPotential()
{}


/// @brief set parameters for distance score between mid points of helices
void
NatbiasHelixPairPotential::set_param_twist( Real w, Real d, Real s )
{
	twist_wts_ = w;
	twist_ = d;
	twist_sigma2_ = s;
}


/// @brief
void
NatbiasHelixPairPotential::score(	SS_Info2_COP const ssinfo, Real & sheet_score ) const
{
	Strands strands( ssinfo->strands() );
	for ( Size istrand=1; istrand<=strands.size(); istrand++ ) {


	}
}  // score


} // ns sspot
} // ns potentials
}	// ns fldsgn
}	// ns protocols
