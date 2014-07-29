// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/screener/RNA_TorsionScreener.cc
/// @brief Screener checking whether the rna torsions are resonable
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

using namespace core;
using namespace core::pose::rna;
using namespace core::chemical::rna;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace screener {

RNA_TorsionScreener::RNA_TorsionScreener() :
	utility::pointer::ReferenceCount(),
	suitename_( new RNA_SuiteName() )
{}

/// @brief screen the pose
bool RNA_TorsionScreener::screen(
	core::pose::Pose const & pose,
	core::Size const suite
) {
	Real const epsilon = angle_conv(
			pose.residue( suite ).mainchain_torsion( EPSILON ) );
	if ( epsilon > suitename_->epsilonmax ) return false;
	if ( epsilon < suitename_->epsilonmin ) return false;

	Real const zeta = angle_conv(
			pose.residue( suite ).mainchain_torsion( ZETA ) );
	if ( zeta > suitename_->zetamax ) return false;
	if ( zeta < suitename_->zetamin ) return false;

	Real const alpha = angle_conv(
			pose.residue( suite + 1 ).mainchain_torsion( ZETA ) );
	if ( alpha > suitename_->alphamax ) return false;
	if ( alpha < suitename_->alphamin ) return false;

	Real const beta = angle_conv(
			pose.residue( suite + 1 ).mainchain_torsion( BETA ) );
	if ( beta > suitename_->betamax ) return false;
	if ( beta < suitename_->betamin ) return false;

	Real const gamma = angle_conv(
			pose.residue( suite + 1 ).mainchain_torsion( GAMMA ) );
	if ( gamma >= suitename_->gammapmin && gamma <= suitename_->gammapmax )
			return true;
	if ( gamma >= suitename_->gammatmin && gamma <= suitename_->gammatmax )
			return true;
	if ( gamma >= suitename_->gammammin && gamma <= suitename_->gammammax )
			return true;

	return false;
}
/////////////////////////////////////////////////////////////////////
Real RNA_TorsionScreener::angle_conv( Real const input ) {
	Real angle =  numeric::principal_angle_degrees( input );
	if ( angle >= 360 ) angle -= 360;
	if ( angle < 0 ) angle += 360;
	return angle;
}
/////////////////////////////////////////////////////////////////////
} //screener
} //sampler
} //stepwise
} //protocols
