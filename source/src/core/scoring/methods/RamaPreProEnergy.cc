// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RamaPreProEnergy.cc
/// @brief  A variation on the Ramachandran scorefunction that has separate probability tables for residues that precede prolines.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- not the original author, but modified this to work with D-amino acids, BACKBONE_AA amino acids, and cyclic
/// geometry.  Refactored greatly to support arbitrary heteropolymer building blocks with any number of mainchain torsions.

// Unit headers
#include <core/scoring/methods/RamaPreProEnergy.hh>
#include <core/scoring/methods/RamaPreProEnergyCreator.hh>
#include <core/scoring/RamaPrePro.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AA.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/PolymerBondedEnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/angle.functions.hh>
#include <core/scoring/DerivVectorPair.hh>

// C++ headers
#include <iostream>
#include <utility/vector1.hh>
#include <core/pose/PDBInfo.hh>

namespace core {
namespace scoring {
namespace methods {


static THREAD_LOCAL basic::Tracer TR( "core.scoring.RamaPreProEnergy" );

//////////////////////
/// EnergyMethod Creator
methods::EnergyMethodOP
RamaPreProEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & /*options*/
) const {
	return methods::EnergyMethodOP( new RamaPreProEnergy( ) );
}

ScoreTypes
RamaPreProEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rama_prepro );
	return sts;
}


RamaPreProEnergy::RamaPreProEnergy( ) :
	parent( methods::EnergyMethodCreatorOP( new RamaPreProEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_RamaPrePro() )
{}

EnergyMethodOP
RamaPreProEnergy::clone() const {
	return EnergyMethodOP( new RamaPreProEnergy( *this ) );
}


void
RamaPreProEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	create_long_range_energy_container( pose, core::scoring::rama_prepro, long_range_type() );
}

bool
RamaPreProEnergy::defines_residue_pair_energy(
	pose::Pose const &pose,
	Size rsd1,
	Size rsd2
) const {
	bool res1_is_lo(false), res2_is_lo(false);
	determine_lo_and_hi_residues( pose, rsd1, rsd2, res1_is_lo, res2_is_lo );
	runtime_assert_string_msg( !(res1_is_lo && res2_is_lo ), "Error in core::scoring::methods::RamaPreProEnergy::defines_residue_pair_energy(): The RamaPrePro term is incompatible with cyclic dipeptides (as is most of the rest of Rosetta)." );

	return (res1_is_lo || res2_is_lo);
}

methods::LongRangeEnergyType
RamaPreProEnergy::long_range_type() const { return methods::ramaprepro_lr; }

void
RamaPreProEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using namespace numeric;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;
	if ( !rsd1.is_bonded(rsd2) ) return;

	bool res1_is_lo(false), res2_is_lo(false);

	determine_lo_and_hi_residues( pose, rsd1.seqpos(), rsd2.seqpos(), res1_is_lo, res2_is_lo );

	if ( !(res1_is_lo || res2_is_lo) ) return;
	runtime_assert_string_msg( !(res1_is_lo && res2_is_lo ), "Error in core::scoring::methods::RamaPreProEnergy::residue_pair_energy(): The RamaPrePro term is incompatible with cyclic dipeptides (as is most of the rest of Rosetta)." );

	conformation::Residue const &res_lo = (res1_is_lo) ? rsd1 : rsd2;
	conformation::Residue const &res_hi = (res2_is_lo) ? rsd1 : rsd2;

	if ( res_lo.has_variant_type(core::chemical::CUTPOINT_LOWER) ) return;
	if ( res_hi.has_variant_type(core::chemical::CUTPOINT_UPPER) ) return;
	if ( res_lo.is_terminus() || !res_lo.has_lower_connect() || !res_lo.has_upper_connect() ) return; // Rama not defined.  Note that is_terminus() checks for UPPER_TERMINUS or LOWER_TERMINUS variant types; it knows nothing about sequence position.

	utility::vector1 < core::Real > mainchain_torsions( res_lo.type().mainchain_atoms().size() - 1 );
	for ( core::Size i=1, imax=mainchain_torsions.size(); i<=imax; ++i ) mainchain_torsions[i] = nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(i) );
	Real rama_score;
	utility::vector1 < core::Real > gradient; //Dummy argument for below.
	potential_.eval_rpp_rama_score( res_lo.type().get_self_ptr(), res_hi.aa(), mainchain_torsions, rama_score, gradient, false /*Don't return gradient*/ );

	emap[ rama_prepro ] += rama_score;
}


Real
RamaPreProEnergy::eval_intraresidue_dof_derivative(
	conformation::Residue const & res_lo,
	ResSingleMinimizationData const & /*min_data*/,
	id::DOF_ID const & /*dof_id*/,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights
) const
{
	using namespace numeric;

	//std::cerr << "!" << std::endl;

	Real deriv(0.0);
	if ( tor_id.valid() && tor_id.type() == id::BB ) {
		if ( tor_id.rsd() != res_lo.seqpos() ) return 0.0;
		if ( res_lo.is_terminus() || !res_lo.has_lower_connect() || !res_lo.has_upper_connect() ) return 0.0;
		if ( pose.fold_tree().is_cutpoint( res_lo.seqpos() ) ) return 0.0;

		conformation::Residue const &res_hi( pose.residue( res_lo.residue_connection_partner( res_lo.upper_connect().index() ) ) );

		if ( tor_id.torsion() == res_lo.type().mainchain_atoms().size() ) return 0.0; //No derivative for omega, the inter-residue torsion.

		utility::vector1 < core::Real > mainchain_torsions( res_lo.type().mainchain_atoms().size() - 1 );
		for ( core::Size i=1, imax=mainchain_torsions.size(); i<=imax; ++i ) mainchain_torsions[i] = nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(i) );
		Real rama_score; //Dummy argument to below.
		utility::vector1 < core::Real > gradient; //Dummy argument for below.
		potential_.eval_rpp_rama_score( res_lo.type().get_self_ptr(), res_hi.aa(), mainchain_torsions, rama_score, gradient, true /*Do return gradient*/ );

		debug_assert( gradient.size() == res_lo.type().mainchain_atoms().size() - 1 );
		deriv = gradient[tor_id.torsion()];
	}

	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return weights[ rama_prepro ] * numeric::conversions::degrees( deriv );
}

core::Size
RamaPreProEnergy::version() const {
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
