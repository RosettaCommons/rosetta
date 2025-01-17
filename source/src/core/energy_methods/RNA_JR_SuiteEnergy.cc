// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_JR_SuiteEnergy.cc
/// @brief  Jane Richardson style RNA Suite Energy
/// @author Fang-Chieh Chou

// Unit Headers
#include <core/energy_methods/RNA_JR_SuiteEnergy.hh>
#include <core/energy_methods/RNA_JR_SuiteEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>

// Project Headers
#include <core/chemical/rna/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

using namespace core::chemical;
using namespace core::chemical::rna;
using namespace core::pose::rna;

/////////////////////////////////////////////////////////////////////////////////////////////////
//
// This makes use of the RNA_SuiteName assignment framework in core/pose/rna/.
//
// It does *not* have derivatives, and the score term available is
//
//   0.0             (for dist7 < 1.0)
//   dist7 - 1.0     (for dist7 > 1.0)
//
// Further development of differential score terms is in RNA_SuiteEnergy.cc &
//  RNA_SuitePotential.cc, including differentiable suiteness. We could probably
//  reconstitute this function in there as well.
//
/////////////////////////////////////////////////////////////////////////////////////////////////


namespace core {
namespace energy_methods {

/// @details This must return a fresh instance,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RNA_JR_SuiteEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const { return utility::pointer::make_shared< RNA_JR_SuiteEnergy >(); }

core::scoring::ScoreTypes
RNA_JR_SuiteEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.emplace_back( rna_jr_suite );
	return sts;
}

/// ctor
RNA_JR_SuiteEnergy::RNA_JR_SuiteEnergy() :
	parent( utility::pointer::make_shared< RNA_JR_SuiteEnergyCreator >() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
RNA_JR_SuiteEnergy::clone() const { return utility::pointer::make_shared< RNA_JR_SuiteEnergy >(); }

///////////////////////////////////////////////////////////////////////////////
void
RNA_JR_SuiteEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	using namespace core::id;
	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd2.has_variant_type( REPLONLY ) ) return;

	Size const rsdnum1( rsd1.seqpos() ), rsdnum2( rsd2.seqpos() );
	if ( rsdnum1 + 1 != rsdnum2 ) return;

	utility::vector1< TorsionID > torsion_ids;
	torsion_ids.emplace_back( rsdnum1, id::BB, DELTA );
	torsion_ids.emplace_back( rsdnum1, id::BB, EPSILON );
	torsion_ids.emplace_back( rsdnum1, id::BB, ZETA );
	torsion_ids.emplace_back( rsdnum2, id::BB, ALPHA );
	torsion_ids.emplace_back( rsdnum2, id::BB, BETA );
	torsion_ids.emplace_back( rsdnum2, id::BB, GAMMA );
	torsion_ids.emplace_back( rsdnum2, id::BB, DELTA );

	utility::fixedsizearray1<Real, 7> torsions;

	for ( Size i = 1; i <= torsion_ids.size(); ++i ) {
		if ( !is_torsion_valid( pose, torsion_ids[i] ) ) return;
		torsions[i] = pose.torsion( torsion_ids[i] );
	}
	RNA_SuiteAssignment assignment( suitename_.assign( torsions ) );
	Real const dist( assignment.dist );
	if ( dist > 1 ) emap[ core::scoring::rna_jr_suite ] += dist - 1;
}

} //scoring
} //core

