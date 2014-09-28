// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/TorsionFeature.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/TorsionFeature.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>

#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace features {

TorsionFeature::TorsionFeature() :
	torsion_bin_( INVALID )
{}

TorsionFeature::TorsionFeature( TorsionFeature const & other ) :
	ResidueFeature(),
	torsion_bin_( other.torsion_bin() )
{}

TorsionFeature::TorsionFeature( TorsionBin bin ) :
	torsion_bin_( bin )
{}

utility::vector1< ResidueFeatureOP >
TorsionFeature::values_from_pose( core::pose::Pose & pose ) const {
	using core::Real;
	using utility::vector1;
	vector1< TorsionFeatureOP > features;

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		vector1< Real > const & torsions( pose.residue(ii).mainchain_torsions() );
		Real phi   = torsions[1];
		Real psi   = torsions[2];
		Real omega = torsions[3];
		features.push_back(
			protocols::comparative_modeling::features::TorsionFeatureOP( new TorsionFeature( torsion2big_bin( phi, psi, omega ) ) )
		);
	}

	return features;
}

std::string TorsionFeature::type() const {
	return "torsion_bin";
}

ResidueFeatureOP TorsionFeature::clone() const {
	TorsionFeatureOP copy( new TorsionFeature( *this ) );
	return copy;
}

TorsionBin TorsionFeature::torsion_bin() const {
	return torsion_bin_;
}

TorsionBin
TorsionFeature::torsion2big_bin(
	core::Real const phi,
	core::Real const psi,
	core::Real const omega
)
{
	if ( std::abs( omega ) < 90.0 ) {
		return O; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100.0 < psi && psi <= 100.0 ) {
			return G; // alpha-L
		} else {
			return E; // E
		}
	} else {
		if ( -125.0 < psi && psi <= 50.0 ) {
			return A; // helical
		} else {
			return B; // extended
		}
	}
	return X;
}

} // features
} // comparative_modeling
} // protocols
