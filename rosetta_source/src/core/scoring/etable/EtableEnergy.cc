// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


// Unit headers
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/EtableEnergyCreator.hh>

// Package headers
#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>

//Auto Headers
#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/trie.functions.hh>


namespace core {
namespace scoring {
namespace etable {


/// @details This must return a fresh instance of the EtableEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
EtableEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new EtableEnergy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
}

ScoreTypes
EtableEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_atr );
	sts.push_back( fa_rep );
	sts.push_back( fa_sol );
	sts.push_back( fa_intra_atr );
	sts.push_back( fa_intra_rep );
	sts.push_back( fa_intra_sol );
	return sts;
}


/// construction with an etable
EtableEnergy::EtableEnergy(
	Etable const & etable_in,
	methods::EnergyMethodOptions const & options
) :
	BaseEtableEnergy< EtableEnergy > ( new EtableEnergyCreator, etable_in, options, fa_atr, fa_rep, fa_sol ),
	using_interres_scoretypes_( true )
{
}


void
EtableEnergy::setup_for_scoring_( pose::Pose const &pose, scoring::ScoreFunction const& ) const
{
	if (pose.total_residue()) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable_.atom_set() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}
}

bool
EtableEnergy::defines_intrares_energy(
	EnergyMap const & weights
) const
{
	return ( weights[ fa_intra_atr ] != 0 ||
		weights[ fa_intra_rep ] != 0 ||
		weights[ fa_intra_sol ] != 0 );
}

/// @brief
void
EtableEnergy::eval_intrares_energy(
	conformation::Residue const & res,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( res.has_variant_type( "REPLONLY" ) ){
			return;
	}

	EnergyMap tbemap;
	if ( pose.energies().use_nblist() ) return; // intraresidue atom pairs present in neighborlist, evaluated during finalize

	if ( using_interres_scoretypes_ ) {
		set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
		using_interres_scoretypes_ = false;
	}
	assert( rep_scoretype() == fa_intra_rep );
	count_pair::CountPairFunctionOP cpfxn = get_intrares_countpair( res, pose, sfxn );
	cpfxn->residue_atom_pair_energy( res, res, *this, tbemap );
	emap[ fa_intra_atr ] = tbemap[ fa_intra_atr ];
	emap[ fa_intra_rep ] = tbemap[ fa_intra_rep ];
	emap[ fa_intra_sol ] = tbemap[ fa_intra_sol ];
}
core::Size
EtableEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace etable
} // namespace scoring
} // namespace core
