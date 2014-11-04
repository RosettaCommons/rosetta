// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/FingerprintMultifunc.hh
/// @brief  Fingerprint multifunction class
/// @author Ragul Gowthaman

/// Unit headers
#include <protocols/pockets/FingerprintMultifunc.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>

#include <cmath>

#include <utility/vector1.hh>


namespace protocols{
namespace pockets {

FingerprintMultifunc::FingerprintMultifunc(
	NonPlaidFingerprint & nfp_in,
	PlaidFingerprint & pfp_in,
	core::Real const & missing_point_weight,
	core::Real const & steric_weight,
	core::Real const & extra_point_weight,
	core::Size const nconformers
) :
	nfp_( nfp_in ),
	pfp_( pfp_in ),
  missing_pt_(missing_point_weight),
	steric_(steric_weight),
  extra_pt_(extra_point_weight),
  nconformers_(nconformers)
{}

core::Real
FingerprintMultifunc::operator ()( core::optimization::Multivec const & vars ) const {
	numeric::xyzVector<core::Real> origin_offset;
	origin_offset.x() = vars[1];
	origin_offset.y() = vars[2];
	origin_offset.z() = vars[3];

	core::Size const conformer=((core::Size)(floor(vars[7])) % nconformers_);
	pfp_.move_ligand_and_update_rhos_( nfp_, origin_offset, vars[4], vars[5], vars[6], conformer );
	core::Real const score = pfp_.fp_compare( nfp_, missing_pt_, steric_, extra_pt_ );
	return score;
}

void
FingerprintMultifunc::dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const
{

	//	std::cout<< "Can't compute derivates of FingerprintMultifunc" << std::endl;
	//	exit(1);

	numeric::xyzVector<core::Real> origin_offset;
	origin_offset.x() = vars[1];
	origin_offset.y() = vars[2];
	origin_offset.z() = vars[3];

	pfp_.move_ligand_and_update_rhos_( nfp_, origin_offset, vars[4], vars[5], vars[6], (core::Size)vars[7] );

	pfp_.fp_compare_deriv( nfp_, missing_pt_, steric_, extra_pt_, dE_dvars[1], dE_dvars[2], dE_dvars[3], dE_dvars[4], dE_dvars[5], dE_dvars[6] );

	return;

}

/// @details Useful debugging code that can be re-enabled by changing the boolean
/// variables at the top of this file.
void
FingerprintMultifunc::dump( core::optimization::Multivec const & vars ) const {
	std::cout<< "In FingerprintMultifunc vars are " << vars[1] << ' ' << vars[2] << ' ' << vars[3] << ' ' << vars[4] << ' ' << vars[5] << ' ' << vars[6] << ' ' << std::endl;
}

} // namespace pockets
} // namespace protocols

