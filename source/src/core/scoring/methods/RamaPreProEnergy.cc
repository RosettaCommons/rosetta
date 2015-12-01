// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RamaPreProEnergy.cc
/// @brief
/// @author

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
	using namespace methods;

	// create LR energy container
	LongRangeEnergyType const & lr_type( long_range_type() );
	Energies & energies( pose.energies() );
	bool create_new_lre_container( false );

	if ( energies.long_range_container( lr_type ) == 0 ) {
		create_new_lre_container = true;
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( lr_type );
		PeptideBondedEnergyContainerOP dec( utility::pointer::static_pointer_cast< core::scoring::PeptideBondedEnergyContainer > ( lrc ) );
		Size nres = pose.total_residue();
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
		if ( dec->size() != nres ) {
			create_new_lre_container = true;
		}
	}

	if ( create_new_lre_container ) {
		Size nres = pose.total_residue();
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			nres = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
		utility::vector1< ScoreType > s_types;
		s_types.push_back( rama_prepro );
		LREnergyContainerOP new_dec( new PeptideBondedEnergyContainer( nres, s_types ) );
		energies.set_long_range_container( lr_type, new_dec );
	}
}

bool
RamaPreProEnergy::defines_residue_pair_energy(
	pose::Pose const &,
	Size res1,
	Size res2
) const {
	return ( res1 == (res2+1) || res1 == (res2-1) );
}

methods::LongRangeEnergyType
RamaPreProEnergy::long_range_type() const { return methods::rama2b_lr; }

void
RamaPreProEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using namespace numeric;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) { return; }
	if ( !rsd1.is_bonded(rsd2) ) { return; }

	conformation::Residue const &res_lo = (rsd1.seqpos()<rsd2.seqpos()) ? rsd1 : rsd2;
	conformation::Residue const &res_hi = (rsd1.seqpos()<rsd2.seqpos()) ? rsd2 : rsd1;

	if ( res_lo.has_variant_type(core::chemical::CUTPOINT_LOWER) ) { return; }
	if ( res_hi.has_variant_type(core::chemical::CUTPOINT_UPPER) ) { return; }
	if ( res_lo.is_terminus() ) { return; } // rama not defined

	Real const phi ( nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(1)));
	Real const psi ( nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(2)));
	Real rama_score, drpp_dphi, drpp_dpsi;
	potential_.eval_rpp_rama_score( res_lo.aa(), res_hi.aa(), phi, psi, rama_score, drpp_dphi, drpp_dpsi );

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
		core::Size i = tor_id.rsd();
		if ( i==1 || i==pose.total_residue() ) return 0.0;

		conformation::Residue const &res_hi( pose.residue( i+1 ) );

		if ( pose.fold_tree().is_cutpoint( res_lo.seqpos() ) ) { return 0.0; }
		if ( res_lo.is_lower_terminus() ) { return 0.0; } // rama not defined

		if ( res_lo.is_protein() && tor_id.torsion() <= 2 && ! res_lo.is_terminus() ) {
			Real const phi ( nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(1)));
			Real const psi ( nonnegative_principal_angle_degrees( res_lo.mainchain_torsion(2)));
			Real rama_score, drpp_dphi, drpp_dpsi;
			potential_.eval_rpp_rama_score( res_lo.aa(), res_hi.aa(), phi, psi, rama_score, drpp_dphi, drpp_dpsi );

			deriv = ( tor_id.torsion() == 1 ? drpp_dphi : drpp_dpsi );
		}
	}
	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ rama_prepro ] * deriv );
}

core::Size
RamaPreProEnergy::version() const {
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
