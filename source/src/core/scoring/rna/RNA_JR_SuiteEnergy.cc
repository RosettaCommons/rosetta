// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_JR_SuiteEnergy.cc
/// @brief  Jane Richardson style RNA Suite Energy
/// @author Fang-Chieh Chou

// Unit Headers
#include <core/scoring/rna/RNA_JR_SuiteEnergy.hh>
#include <core/scoring/rna/RNA_JR_SuiteEnergyCreator.hh>

// Package Headers
#include <core/scoring/ScoringManager.hh>
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

using namespace core::chemical::rna;
using namespace core::pose::rna;

namespace core {
namespace scoring {
namespace rna {
/// @details This must return a fresh instance,
/// never an instance already in use
methods::EnergyMethodOP
RNA_JR_SuiteEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {	return new RNA_JR_SuiteEnergy; }

ScoreTypes
RNA_JR_SuiteEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_jr_suite );
	return sts;
}

/// ctor
RNA_JR_SuiteEnergy::RNA_JR_SuiteEnergy() :
	parent( new RNA_JR_SuiteEnergyCreator )
{}

/// clone
methods::EnergyMethodOP
RNA_JR_SuiteEnergy::clone() const { return new RNA_JR_SuiteEnergy; }

///////////////////////////////////////////////////////////////////////////////
void
RNA_JR_SuiteEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using namespace core::id;

	Size const rsdnum1( rsd1.seqpos() ), rsdnum2( rsd2.seqpos() );
	if ( rsdnum1 + 1 != rsdnum2 ) return;

	utility::vector1<TorsionID> torsion_ids;
	torsion_ids.push_back( TorsionID( rsdnum1, id::BB, DELTA ) );
	torsion_ids.push_back( TorsionID( rsdnum1, id::BB, EPSILON ) );
	torsion_ids.push_back( TorsionID( rsdnum1, id::BB, ZETA ) );
	torsion_ids.push_back( TorsionID( rsdnum2, id::BB, ALPHA ) );
	torsion_ids.push_back( TorsionID( rsdnum2, id::BB, BETA ) );
	torsion_ids.push_back( TorsionID( rsdnum2, id::BB, GAMMA ) );
	torsion_ids.push_back( TorsionID( rsdnum2, id::BB, DELTA ) );

	utility::vector1<Real> torsions;

	for ( Size i = 1; i <= torsion_ids.size(); ++i ) {
		if ( !is_torsion_valid( pose, torsion_ids[i] ) ) return;
		Real const tor( pose.torsion( torsion_ids[i] ) );
		torsions.push_back( tor );
	}
	RNA_SuiteAssignment assignment( suitename_.assign( torsions ) );
	Real const dist( assignment.dist );
	if ( dist > 1 ) emap[ rna_jr_suite ] += dist - 1;
}

} //rna
} //scoring
} //core

