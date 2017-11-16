// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_FullAtomVDW_BasePhosphate
/// @brief  RNA_FullAtomVDW_BasePhosphate energy method class implementation
/// @author Parin Sripakdeevong

// Unit Headers
#include <core/scoring/rna/RNA_FullAtomVDW_BasePhosphate.hh>
#include <core/scoring/rna/RNA_FullAtomVDW_BasePhosphateCreator.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>

// Package Headers
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/types.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/conversions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/NeighborList.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>
// C++

static basic::Tracer TR( "core.scoring.rna.RNA_FullAtomVDW_BasePhosphateEnergy", basic::t_info );

///////////////////////////////////////////////////////////////////////////////////////////
//
// Totally hacky score term, which should be deprecated soon, maybe 2015.
//  See note below from Parin, 2012 -- apparently should not have gone
//  into trunk at all, but was used for CS-ROSETTA-RNA paper (Nature Methods, 2014).
//
// Was in use to handle syn-G base-phosphate hydrogen bonds in RNA.
//
// Should be replaceable by standard fa_atr, fa_rep with
//   PUT_INTRA_INTO_TOTAL tag in .wts files ( or
//   use separate fa_atr_intra_crossover4, etc. terms.)
//
//
// -- rhiju, 2014
//
// AMW: Since not in standard use, not bothering to special-case REPLONLY.
////////////////////////////////////////////////////////////////////////////////////////////

using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace rna {

/// @details This must return a fresh instance of the RNA_FullAtomVDW_BasePhosphateCreator class
/// never an instance already in use
methods::EnergyMethodOP
RNA_FullAtomVDW_BasePhosphateCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new RNA_FullAtomVDW_BasePhosphate( options ) );
}

ScoreTypes
RNA_FullAtomVDW_BasePhosphateCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_intra_RNA_base_phos_atr );
	sts.push_back( fa_intra_RNA_base_phos_rep );
	sts.push_back( fa_intra_RNA_base_phos_sol );
	return sts;
}


/// constructor
RNA_FullAtomVDW_BasePhosphate::RNA_FullAtomVDW_BasePhosphate(
	methods::EnergyMethodOptions const & options
):
	parent( methods::EnergyMethodCreatorOP( new RNA_FullAtomVDW_BasePhosphateCreator ) ),
	options_( options )
{
	etable::Etable const & etable = *(ScoringManager::get_instance()->etable( options ).lock());
	if ( options.analytic_etable_evaluation() ) {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::AnalyticEtableEvaluator( etable ) );
	} else {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::TableLookupEvaluator( etable ) );
	}

	etable_evaluator_->set_scoretypes(
		fa_intra_RNA_base_phos_atr,
		fa_intra_RNA_base_phos_rep,
		fa_intra_RNA_base_phos_sol );
}


RNA_FullAtomVDW_BasePhosphate::~RNA_FullAtomVDW_BasePhosphate() {}


/// clone
methods::EnergyMethodOP
RNA_FullAtomVDW_BasePhosphate::clone() const
{

	return methods::EnergyMethodOP( new RNA_FullAtomVDW_BasePhosphate( *this ) );
}

void
RNA_FullAtomVDW_BasePhosphate::residue_fast_pair_energy_attached_H(
	conformation::Residue const & res1,
	int const atomno1,
	conformation::Residue const & res2,
	Size const atomno2,
	Size const at1hbegin, //at1hbegin and at1hend define a range of hydrogen atom indices -- those h's bound to at1
	Size const at1hend,
	Size const at2hbegin,
	Size const at2hend,
	EnergyMap & emap
) const {
	using conformation::Atom;

	Weight weight( 1.0 );

	Atom const & atom1( res1.atom( atomno1 ) );
	Atom const & atom2( res2.atom( atomno2 ) );

	// Heavy Atom in res1 to Hs in res2
	for ( Size i = at2hbegin; i <= at2hend; ++i ) {
		Atom const & H2( res2.atom( i ) );
		weight = 1.0;
		etable_evaluator_->pair_energy_H_v( atom1, H2, weight, emap );
	}

	// Hs in res1 to heavy Atom and Hs in res2
	for ( Size i = at1hbegin; i <= at1hend; ++i ) {
		Atom const & H1( res1.atom( i ) );
		weight = 1.0;
		// H in res1 to heavy Atom in res2
		etable_evaluator_->pair_energy_H_v( H1, atom2, weight, emap );

		// H in res1 to Hs in res2
		for ( Size j = at2hbegin; j <= at2hend; ++j ) {
			Atom const & H2( res2.atom( j ) );
			weight = 1.0f;
			etable_evaluator_->pair_energy_H_v( H1, H2, weight, emap );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void
RNA_FullAtomVDW_BasePhosphate::residue_energy(
	conformation::Residue const & rsd,
	EnergyMap & emap  ) const {
	using conformation::Atom;

	if ( rsd.is_RNA() == false ) return;

	DistanceSquared dsq;

	//Weight const weight=1.0;
	Real const weight = 1.0;

	// get hydrogen interaction cutoff
	Real const Hydrogen_interaction_cutoff2 = ( etable_evaluator_->hydrogen_interaction_cutoff2() );

	typedef utility::vector1< Size > const & vect;

	vect rhbegin( rsd.attached_H_begin() );
	vect rhend(   rsd.attached_H_end()   );

	Size const rsdnheavyatoms = rsd.nheavyatoms();

	// Atom pairs
	for ( Size i = 1; i <= rsdnheavyatoms; ++i ) {
		Atom const & atom1( rsd.atom( i ) );
		for ( Size j = i + 1; j <= rsdnheavyatoms; ++j ) {
			Atom const & atom2( rsd.atom( j ) );

			if ( !core::chemical::rna::is_base_phosphate_atom_pair( rsd, rsd, i, j ) ) continue;

			etable_evaluator_->atom_pair_energy_v( atom1, atom2, weight, emap, dsq );

			if ( dsq < Hydrogen_interaction_cutoff2 ) {
				residue_fast_pair_energy_attached_H( rsd, i, rsd, j, rhbegin[ i ], rhend[ i ], rhbegin[ j ], rhend[ j ], emap );
			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomVDW_BasePhosphate::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap  ) const {

	return residue_energy( rsd, emap );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_FullAtomVDW_BasePhosphate::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & /*sfxn*/, // needed for non-nblist minimization
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	if ( weights[ fa_intra_RNA_base_phos_sol] != 0.0 ) {
		//Please refer to paragraph at the end of RNA_FullAtomVDW_BasePhosphate::residue_energy for explanation.
		//Again, if you want to implement this term, please ensure your implementation works properly (i.e. perform numerical_derivative_check() and etc),
		//before committing the code to TRUNK!)
		// Parin S. (sripakpa@stanford.edu). Jan 11, 2012
		utility_exit_with_message( "weights[ fa_intra_RNA_base_phos_sol ] != 0.0, but this term is not yet implemented!" );
	}

	Size const seq_num = id.rsd();

	if ( pose.residue( seq_num ).is_RNA() == false ) return;

	conformation::Residue const & rsd = pose.residue( seq_num );

	conformation::Atom const & atom1( rsd.atom( id.atomno() ) );

	Real const cp_weight = 1.0;

	Vector f1, f2;
	for ( Size nbr_atomno = 1; nbr_atomno <= rsd.natoms(); nbr_atomno++ ) {

		if ( rsd.is_virtual( id.atomno() ) ) continue; //Is this necessary?
		if ( rsd.is_virtual( nbr_atomno  ) ) continue; //Is this necessary?
		if ( !core::chemical::rna::is_base_phosphate_atom_pair( rsd, rsd, id.atomno(), nbr_atomno ) ) continue;

		if ( rsd.path_distance( nbr_atomno, id.atomno() ) < 4 ) utility_exit_with_message( "rsd.path_distance( nbr_atomno, id.atomno() ) < 4" );

		conformation::Atom const & atom2( rsd.atom( nbr_atomno ) );

		Real const dE_dR_over_r = etable_evaluator_->eval_dE_dR_over_r_v( atom1, atom2, weights, f1, f2 );

		if ( dE_dR_over_r != 0.0 ) {
			F1 += dE_dR_over_r * cp_weight * f1;
			F2 += dE_dR_over_r * cp_weight * f2;
		}
	}
}

void RNA_FullAtomVDW_BasePhosphate::indicate_required_context_graphs( utility::vector1< bool > & ) const{}

core::Size
RNA_FullAtomVDW_BasePhosphate::version() const
{
	return 2; // After bug fix to actually store energies, rhiju Oct. 2014.
	// return 1; // First version, created by Parin Sripakdeevong (sripakpa@stanford.edu), Jan 2012.
}


} //rna
} //scoring
} //core

