// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ReferenceEnergy.hh
/// @brief  Reference energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/methods/ReferenceEnergy.hh>
#include <core/scoring/methods/ReferenceEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the ReferenceEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ReferenceEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	if ( options.has_method_weights( ref ) ) {
		return methods::EnergyMethodOP( new ReferenceEnergy( options.method_weights( ref ) ) );
	} else {
		return methods::EnergyMethodOP( new ReferenceEnergy );
	}
}

ScoreTypes
ReferenceEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ref );
	return sts;
}


ReferenceEnergy::ReferenceEnergy() :
	parent( methods::EnergyMethodCreatorOP( new ReferenceEnergyCreator ) )
{}

ReferenceEnergy::ReferenceEnergy( utility::vector1< Real > const & aa_weights_in ):
	parent( methods::EnergyMethodCreatorOP( new ReferenceEnergyCreator ) ),
	aa_weights_( aa_weights_in )
{}

ReferenceEnergy::~ReferenceEnergy() {}

EnergyMethodOP
ReferenceEnergy::clone() const
{
	return EnergyMethodOP( new ReferenceEnergy( aa_weights_ ) );
}


/// This is a terrible terrible terrible hack that will do for now.
void
ReferenceEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	using namespace chemical;
	if ( !aa_weights_.empty() ) {
		///
		AA const & aa( rsd.aa() );
		AA const aa2 ( is_d_aminoacid(aa) ? get_l_equivalent(aa) : aa);
		if ( Size(aa2) > aa_weights_.size() ) return;
		emap[ ref ] += aa_weights_[ aa2 ];
		//   if ( rsd.is_DNA() ) {
		//    std::cout << "using dna refE " << aa_weights_[aa] << std::endl;
		//   }
		return;
	}

	// reference weights for RNA...
	if ( rsd.is_RNA() && ( aa_weights_.size() > num_canonical_aas ) ) {
		emap[ ref ] += aa_weights_[ rsd.aa() ];
		return;
	}

	/// else -- use the default reference weights from r++
	if ( rsd.type().aa() > num_canonical_aas ) return;

	switch ( rsd.type().aa() ) {
	case aa_ala : emap[ ref ] +=  0.16; break;
	case aa_cys : emap[ ref ] +=  1.70; break;
	case aa_asp :
		emap[ ref ] += (rsd.type().has_variant_type( chemical::PROTONATED )) ? -0.262 : -0.67;
		break;
	case aa_glu :
		emap[ ref ] += (rsd.type().has_variant_type( chemical::PROTONATED )) ? -0.81 : -0.81;
		break;
	case aa_phe : emap[ ref ] +=  0.63; break;
	case aa_gly : emap[ ref ] +=  -0.17; break;
	case aa_his :
		emap[ ref ] += (rsd.type().has_variant_type( chemical::PROTONATED )) ? 0.288 : 0.56;
		break;
	case aa_ile : emap[ ref ] +=  0.24; break;
	case aa_lys :
		emap[ ref ] += (rsd.type().has_variant_type( chemical::DEPROTONATED )) ? -0.65 : -0.65;
		break;
	case aa_leu : emap[ ref ] +=  -0.10; break;
	case aa_met : emap[ ref ] +=  -0.34; break;
	case aa_asn : emap[ ref ] +=  -0.89; break;
	case aa_pro : emap[ ref ] +=  0.02; break;
	case aa_gln : emap[ ref ] +=  -0.97; break;
	case aa_arg : emap[ ref ] +=  -0.98; break;
	case aa_ser : emap[ ref ] +=  -0.37; break;
	case aa_thr : emap[ ref ] +=  -0.27; break;
	case aa_val : emap[ ref ] +=  0.29; break;
	case aa_trp : emap[ ref ] +=  0.91; break;
	case aa_tyr :
		emap[ ref ] += (rsd.type().has_variant_type( chemical::DEPROTONATED )) ? 0.238 : 0.51;
		break;
	default :
		utility_exit();
		break;
	}

}

///////////////////////////////////////////////////////////////////////////////

/// @brief Returns true if passed a core::chemical::AA corresponding to a
/// D-amino acid, and false otherwise.
bool
ReferenceEnergy::is_d_aminoacid(
	core::chemical::AA const res_aa
) const {
	using namespace core::chemical;
	if ( res_aa >= aa_dal && res_aa <= aa_dty ) return true;
	return false;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief When passed a d-amino acid, returns the l-equivalent.  Returns
/// aa_unk otherwise.
core::chemical::AA
ReferenceEnergy::get_l_equivalent(
	core::chemical::AA const d_aa
) const {
	using namespace core::chemical;
	if ( d_aa==aa_dal ) return aa_ala;
	else if ( d_aa==aa_dcs ) return aa_cys;
	else if ( d_aa==aa_das ) return aa_asp;
	else if ( d_aa==aa_dgu ) return aa_glu;
	else if ( d_aa==aa_dph ) return aa_phe;
	else if ( d_aa==aa_dhi ) return aa_his;
	else if ( d_aa==aa_dil ) return aa_ile;
	else if ( d_aa==aa_dly ) return aa_lys;
	else if ( d_aa==aa_dle ) return aa_leu;
	else if ( d_aa==aa_dme ) return aa_met;
	else if ( d_aa==aa_dan ) return aa_asn;
	else if ( d_aa==aa_dpr ) return aa_pro;
	else if ( d_aa==aa_dgn ) return aa_gln;
	else if ( d_aa==aa_dar ) return aa_arg;
	else if ( d_aa==aa_dse ) return aa_ser;
	else if ( d_aa==aa_dth ) return aa_thr;
	else if ( d_aa==aa_dva ) return aa_val;
	else if ( d_aa==aa_dtr ) return aa_trp;
	else if ( d_aa==aa_dty ) return aa_tyr;

	return aa_unk;
}


Real
ReferenceEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

/// @brief ReferenceEnergy is context independent; indicates that no
/// context graphs are required
void
ReferenceEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
ReferenceEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

