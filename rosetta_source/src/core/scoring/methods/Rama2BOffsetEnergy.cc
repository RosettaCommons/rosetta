// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Rama2BOffsetEnergy.cc
/// @brief
/// @author

// Unit headers
#include <core/scoring/methods/Rama2BOffsetEnergy.hh>
#include <core/scoring/methods/Rama2BOffsetEnergyCreator.hh>
#include <core/scoring/Rama2BOffset.hh>

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

#include <core/scoring/PeptideBondedEnergyContainer.hh>
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

// options
//#include <basic/options/option.hh>
//#include <basic/options/keys/score.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <utility/vector1.hh>
#include <core/pose/PDBInfo.hh>

namespace core {
namespace scoring {
namespace methods {


static basic::Tracer TR("core.scoring.Rama2BOffsetEnergy");

//////////////////////
/// EnergyMethod Creator
methods::EnergyMethodOP
Rama2BOffsetEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & /*options*/
) const {
	return new Rama2BOffsetEnergy( );
}

ScoreTypes
Rama2BOffsetEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rama2b_offset );
	sts.push_back( omega2b_offset );
	sts.push_back( p_aa_pp_offset );
	return sts;
}


Rama2BOffsetEnergy::Rama2BOffsetEnergy( ) :
	parent( new Rama2BOffsetEnergyCreator ),
	potential_( ScoringManager::get_instance()->get_Rama2BOffset() )
{}

Rama2BOffsetEnergy::~Rama2BOffsetEnergy() {}

EnergyMethodOP
Rama2BOffsetEnergy::clone() const {
	return new Rama2BOffsetEnergy( *this );
}


void
Rama2BOffsetEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	using namespace methods;

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		PeptideBondedEnergyContainerOP dec( static_cast< PeptideBondedEnergyContainer * > ( lrc.get() ) );
		Size nres = pose.total_residue();
		if( core::pose::symmetry::is_symmetric(pose) )
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		if ( dec->size() != nres ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		Size nres = pose.total_residue();
		if( core::pose::symmetry::is_symmetric(pose) )
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		utility::vector1< ScoreType > s_types;
		s_types.push_back( rama2b_offset );
		s_types.push_back( omega2b_offset );
		s_types.push_back( p_aa_pp_offset );
		LREnergyContainerOP new_dec = new PeptideBondedEnergyContainer( nres, s_types );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

bool
Rama2BOffsetEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size res1,
	Size res2
) const {
	return ( res1 == (res2+1) || res1 == (res2-1) );
}

methods::LongRangeEnergyType
Rama2BOffsetEnergy::long_range_type() const { return methods::rama2b_lr; }

void
Rama2BOffsetEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using namespace numeric;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) { return; }
	if ( !rsd1.is_bonded(rsd2) ) { return; }

	conformation::Residue const &res_lo = (rsd1.seqpos()<rsd2.seqpos()) ? rsd1 : rsd2;
	conformation::Residue const &res_hi = (rsd1.seqpos()<rsd2.seqpos()) ? rsd2 : rsd1;

	if ( pose.fold_tree().is_cutpoint( res_lo.seqpos() ) ) { return; }

	// hardcoded atom indices assumes N=1, CA=2, C=3
	core::Real psi1 = numeric::dihedral_degrees(
		res_lo.atom(1).xyz(), res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz() );
	core::Real omega2 = numeric::dihedral_degrees(
		res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz() );
	core::Real phi2 = numeric::dihedral_degrees(
		res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), res_hi.atom(3).xyz() );

	//if (psi1<0)   psi1+=360;
	//if (omega2<0) omega2+=360;
	//if (phi2<0)   phi2+=360;

	Real rama_score, omega_score, paapp_score, drb2o_dpsi1, drb2o_domega2, drb2o_dphi2;
	potential_.eval_r2bo_rama_score( res_lo.aa(), res_hi.aa(), psi1, omega2, phi2, rama_score, drb2o_dpsi1, drb2o_domega2, drb2o_dphi2 );
	potential_.eval_r2bo_omega_score( res_lo.aa(), res_hi.aa(), psi1, omega2, phi2, omega_score, drb2o_dpsi1, drb2o_domega2, drb2o_dphi2 );
	potential_.eval_p_aa_pp_score( res_lo.aa(), res_hi.aa(), psi1, phi2, paapp_score, drb2o_dpsi1, drb2o_dphi2 );

	emap[ rama2b_offset ] += rama_score;
	emap[ omega2b_offset ] += omega_score;
	emap[ p_aa_pp_offset ] += paapp_score;
}


void
Rama2BOffsetEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	using namespace numeric;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) { return; }
	if ( !rsd1.is_bonded(rsd2) ) { return; }

	conformation::Residue const &res_lo = (rsd1.seqpos()<rsd2.seqpos()) ? rsd1 : rsd2;
	conformation::Residue const &res_hi = (rsd1.seqpos()<rsd2.seqpos()) ? rsd2 : rsd1;
	utility::vector1< DerivVectorPair > &res_lo_derivs = (rsd1.seqpos()<rsd2.seqpos()) ? r1_atom_derivs : r2_atom_derivs;
	utility::vector1< DerivVectorPair > &res_hi_derivs = (rsd1.seqpos()<rsd2.seqpos()) ? r2_atom_derivs : r1_atom_derivs;

	if ( pose.fold_tree().is_cutpoint( res_lo.seqpos() ) ) { return; }

	// hardcoded atom indices assumes N=1, CA=2, C=3
	core::Real psi1 = numeric::dihedral_degrees(
		res_lo.atom(1).xyz(), res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz() );
	core::Real omega2 = numeric::dihedral_degrees(
		res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz() );
	core::Real phi2 = numeric::dihedral_degrees(
		res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), res_hi.atom(3).xyz() );

	//if (psi1<0)   psi1+=360;
	//if (omega2<0) omega2+=360;
	//if (phi2<0)   phi2+=360;

	Real rama_score, omega_score, paapp_score;
	Real drama_dpsi1, drama_domega2, drama_dphi2, domega_dpsi1, domega_domega2, domega_dphi2, dpaapp_dpsi1, dpaapp_dphi2;
	potential_.eval_r2bo_rama_score ( res_lo.aa(), res_hi.aa(), psi1, omega2, phi2, rama_score, drama_dpsi1, drama_domega2, drama_dphi2 );
	potential_.eval_r2bo_omega_score( res_lo.aa(), res_hi.aa(), psi1, omega2, phi2, omega_score, domega_dpsi1, domega_domega2, domega_dphi2 );
	potential_.eval_p_aa_pp_score( res_lo.aa(), res_hi.aa(), psi1, phi2, paapp_score, dpaapp_dpsi1, dpaapp_dphi2 );

	Real drb2o_dpsi1, drb2o_domega2, drb2o_dphi2;
	drb2o_dpsi1  = weights[ rama2b_offset ] * numeric::conversions::degrees(drama_dpsi1);
	drb2o_dpsi1 += weights[ omega2b_offset ] * numeric::conversions::degrees(domega_dpsi1);
	drb2o_dpsi1 += weights[ p_aa_pp_offset ] * numeric::conversions::degrees(dpaapp_dpsi1);
	drb2o_domega2  = weights[ rama2b_offset ] * numeric::conversions::degrees(drama_domega2);
	drb2o_domega2 += weights[ omega2b_offset ] * numeric::conversions::degrees(domega_domega2);
	drb2o_dphi2  = weights[ rama2b_offset ] * numeric::conversions::degrees(drama_dphi2);
	drb2o_dphi2 += weights[ omega2b_offset ] * numeric::conversions::degrees(domega_dphi2);
	drb2o_dphi2 += weights[ p_aa_pp_offset ] * numeric::conversions::degrees(dpaapp_dphi2);

	// torsion->cart
	Vector f1, f2;
	Real phi;
	{
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_lo.atom(1).xyz(), res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), phi, f1, f2 );
			res_lo_derivs[ 1 ].f1() += drb2o_dpsi1 * f1;
			res_lo_derivs[ 1 ].f2() += drb2o_dpsi1 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_lo.atom(1).xyz(), res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), phi, f1, f2 );
			res_lo_derivs[ 2 ].f1() += drb2o_dpsi1 * f1;
			res_lo_derivs[ 2 ].f2() += drb2o_dpsi1 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), res_lo.atom(2).xyz(), res_lo.atom(1).xyz(), phi, f1, f2 );
			res_lo_derivs[ 3 ].f1() += drb2o_dpsi1 * f1;
			res_lo_derivs[ 3 ].f2() += drb2o_dpsi1 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), res_lo.atom(2).xyz(), res_lo.atom(1).xyz(), phi, f1, f2 );
			res_hi_derivs[ 1 ].f1() += drb2o_dpsi1 * f1;
			res_hi_derivs[ 1 ].f2() += drb2o_dpsi1 * f2;
	}
	{
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), phi, f1, f2 );
			res_lo_derivs[ 2 ].f1() += drb2o_domega2 * f1;
			res_lo_derivs[ 2 ].f2() += drb2o_domega2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_lo.atom(2).xyz(), res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), phi, f1, f2 );
			res_lo_derivs[ 3 ].f1() += drb2o_domega2 * f1;
			res_lo_derivs[ 3 ].f2() += drb2o_domega2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_hi.atom(2).xyz(), res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), res_lo.atom(2).xyz(), phi, f1, f2 );
			res_hi_derivs[ 1 ].f1() += drb2o_domega2 * f1;
			res_hi_derivs[ 1 ].f2() += drb2o_domega2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_hi.atom(2).xyz(), res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), res_lo.atom(2).xyz(), phi, f1, f2 );
			res_hi_derivs[ 2 ].f1() += drb2o_domega2 * f1;
			res_hi_derivs[ 2 ].f2() += drb2o_domega2 * f2;
	}
	{
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), res_hi.atom(3).xyz(), phi, f1, f2 );
			res_lo_derivs[ 3 ].f1() += drb2o_dphi2 * f1;
			res_lo_derivs[ 3 ].f2() += drb2o_dphi2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_lo.atom(3).xyz(), res_hi.atom(1).xyz(), res_hi.atom(2).xyz(), res_hi.atom(3).xyz(), phi, f1, f2 );
			res_hi_derivs[ 1 ].f1() += drb2o_dphi2 * f1;
			res_hi_derivs[ 1 ].f2() += drb2o_dphi2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv(
				res_hi.atom(3).xyz(), res_hi.atom(2).xyz(), res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), phi, f1, f2 );
			res_hi_derivs[ 2 ].f1() += drb2o_dphi2 * f1;
			res_hi_derivs[ 2 ].f2() += drb2o_dphi2 * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv(
				res_hi.atom(3).xyz(), res_hi.atom(2).xyz(), res_hi.atom(1).xyz(), res_lo.atom(3).xyz(), phi, f1, f2 );
			res_hi_derivs[ 3 ].f1() += drb2o_dphi2 * f1;
			res_hi_derivs[ 3 ].f2() += drb2o_dphi2 * f2;
	}
}

core::Size
Rama2BOffsetEnergy::version() const {
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
