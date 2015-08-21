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
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <sstream>


namespace core {
namespace scoring {
namespace etable {


/// @details This must return a fresh instance of the EtableEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
EtableEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	/// The command line option needs to override the EnergyMethodOptions because an Etable initialized with
	/// the analytic_etable_evaluation flag on will not have allocated the large etables necessary for the
	/// TableLookupEtableEnergy class.
	// NEEDS FIXING?: using reference of temporary locked weak_ptr
	if ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] || options.analytic_etable_evaluation() ) {
		return methods::EnergyMethodOP( new AnalyticEtableEnergy( *( ScoringManager::get_instance()->etable( options ).lock() ), options, false /*do_classic_intrares*/ ) );
	} else {
		return methods::EnergyMethodOP( new TableLookupEtableEnergy( *( ScoringManager::get_instance()->etable( options ).lock() ), options, false /*do_classic_intrares*/ ) );
	}
}

ScoreTypes
EtableEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_atr );
	sts.push_back( fa_rep );
	sts.push_back( fa_sol );
	// xover4 terms calculate intra for atoms with path distance >= 4 (i.e., not involved in same dihedral). -- rhiju
	sts.push_back( fa_intra_atr_xover4 );
	sts.push_back( fa_intra_rep_xover4 );
	sts.push_back( fa_intra_sol_xover4 );
	//  note that 'classic' fa_intra_atr, fa_intra_rep, and fa_intra_sol terms are now in EtableClassicIntraEnergy. -- rhiju
	return sts;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
methods::EnergyMethodOP
EtableClassicIntraEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	if ( basic::options::option[ basic::options::OptionKeys::score::analytic_etable_evaluation ] || options.analytic_etable_evaluation() ) {
		return methods::EnergyMethodOP( new AnalyticEtableEnergy( *( ScoringManager::get_instance()->etable( options ).lock() ), options, true /*do_classic_intrares*/ ) );
	} else {
		return methods::EnergyMethodOP( new TableLookupEtableEnergy( *( ScoringManager::get_instance()->etable( options ).lock() ), options, true /*do_classic_intrares*/ ) );
	}
}

ScoreTypes
EtableClassicIntraEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_atr_dummy );
	sts.push_back( fa_rep_dummy );
	sts.push_back( fa_sol_dummy );
	sts.push_back( fa_intra_atr );
	sts.push_back( fa_intra_rep );
	sts.push_back( fa_intra_sol );
	return sts;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

/// construction with an etable
TableLookupEtableEnergy::TableLookupEtableEnergy(
	Etable const & etable_in,
	methods::EnergyMethodOptions const & options,
	bool const do_classic_intrares
)
:
	BaseEtableEnergy< TableLookupEtableEnergy > ( do_classic_intrares ?
	methods::EnergyMethodCreatorOP( new EtableClassicIntraEnergyCreator ) :
	methods::EnergyMethodCreatorOP( new EtableEnergyCreator ),
	etable_in, options, do_classic_intrares ),
	intrares_evaluator_( etable_in ),
	interres_evaluator_( etable_in )
{
	if ( do_classic_intrares ) {
		interres_evaluator_.set_scoretypes( fa_atr_dummy, fa_rep_dummy, fa_sol_dummy );
		//  interres_evaluator_.set_scoretypes( fa_atr, fa_rep, fa_sol ); // will get zeros.
		intrares_evaluator_.set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
	} else {
		interres_evaluator_.set_scoretypes( fa_atr, fa_rep, fa_sol );
		if ( put_intra_into_total() ) {
			intrares_evaluator_.set_scoretypes( fa_atr, fa_rep, fa_sol );
		} else {
			intrares_evaluator_.set_scoretypes( fa_intra_atr_xover4, fa_intra_rep_xover4, fa_intra_sol_xover4 );
		}
	}
}

TableLookupEtableEnergy::TableLookupEtableEnergy(
	TableLookupEtableEnergy const & src
) :
	BaseEtableEnergy< TableLookupEtableEnergy >( src ),
	intrares_evaluator_( src.intrares_evaluator_ ),
	interres_evaluator_( src.interres_evaluator_ )
{
}

methods::EnergyMethodOP
TableLookupEtableEnergy::clone() const {
	return methods::EnergyMethodOP( new TableLookupEtableEnergy( *this ) );
}

void
TableLookupEtableEnergy::setup_for_scoring_( pose::Pose const &pose, scoring::ScoreFunction const& ) const
{
	if ( pose.total_residue() ) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable().atom_set().lock() ) {
			std::stringstream err_msg;
			err_msg
				<< "Illegal attempt to score with non-identical atom set between pose and etable" << std::endl
				<< "\tpose   atom_type_set: '" << pose.residue(1).type().atom_type_set_ptr()->name() << "'" << std::endl
				<< "\tetable atom_type_set: '" << etable().atom_set().lock()->name() << "'" << std::endl;
			utility_exit_with_message( err_msg.str());
		}
	}
}


bool
TableLookupEtableEnergy::defines_intrares_energy(
	EnergyMap const & weights
) const
{
	if ( do_classic_intrares() ) {
		return ( weights[ fa_intra_atr ] != 0 ||
			weights[ fa_intra_rep ] != 0 ||
			weights[ fa_intra_sol ] != 0 );
	}
	if ( put_intra_into_total() ) {
		return ( weights[ fa_atr ] != 0 ||
			weights[ fa_rep ] != 0 ||
			weights[ fa_sol ] != 0 );
	}
	return ( weights[ fa_intra_atr_xover4 ] != 0 ||
		weights[ fa_intra_rep_xover4 ] != 0 ||
		weights[ fa_intra_sol_xover4 ] != 0 );
}

/// @brief
void
TableLookupEtableEnergy::eval_intrares_energy(
	conformation::Residue const & res,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( res.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	EnergyMap tbemap;
	if ( pose.energies().use_nblist() ) return; // intraresidue atom pairs present in neighborlist, evaluated during finalize

	count_pair::CountPairFunctionOP cpfxn = get_intrares_countpair( res, pose, sfxn );
	//cpfxn->residue_atom_pair_energy( res, res, *this, tbemap );
	intrares_evaluator_.residue_atom_pair_energy( res, res, *cpfxn, tbemap );
	emap[ intrares_evaluator_.st_atr() ] = tbemap[ intrares_evaluator_.st_atr() ];
	emap[ intrares_evaluator_.st_rep() ] = tbemap[ intrares_evaluator_.st_rep() ];
	emap[ intrares_evaluator_.st_sol() ] = tbemap[ intrares_evaluator_.st_sol() ];

}

/// @details
/// Version 2: apl - 2012/06/14 -- Etable smoothing change for LK sol.  Placing the
///   start point for the spline that ramps the LJatr term to zero at the minimum of
///   a) the VDW sum and b) (max_dis - 1.5), which, for famaxdis of 6 is 4.5 A.  This
///   change effects the Br/Br and Br/I energies but no other atom-type pairs.
core::Size
TableLookupEtableEnergy::version() const
{
	return 2;
	//return 1; // Initial versioning
}


/////////////////////// Analytic Etable Energy Evaluation /////////////////////////////////


/// construction with an etable
AnalyticEtableEnergy::AnalyticEtableEnergy(
	Etable const & etable_in,
	methods::EnergyMethodOptions const & options,
	bool const do_classic_intrares
) :
	BaseEtableEnergy< AnalyticEtableEnergy > ( do_classic_intrares ?
	methods::EnergyMethodCreatorOP( new EtableClassicIntraEnergyCreator ) :
	methods::EnergyMethodCreatorOP( new EtableEnergyCreator ),
	etable_in, options, do_classic_intrares ),
	intrares_evaluator_( etable_in ),
	interres_evaluator_( etable_in )
{
	if ( do_classic_intrares ) {
		interres_evaluator_.set_scoretypes( fa_atr_dummy, fa_rep_dummy, fa_sol_dummy );
		intrares_evaluator_.set_scoretypes( fa_intra_atr, fa_intra_rep, fa_intra_sol );
	} else {
		interres_evaluator_.set_scoretypes( fa_atr, fa_rep, fa_sol );
		if ( put_intra_into_total() ) {
			intrares_evaluator_.set_scoretypes( fa_atr, fa_rep, fa_sol );
		} else {
			intrares_evaluator_.set_scoretypes( fa_intra_atr_xover4, fa_intra_rep_xover4, fa_intra_sol_xover4 );
		}
	}
}

AnalyticEtableEnergy::AnalyticEtableEnergy(
	AnalyticEtableEnergy const & src
) :
	BaseEtableEnergy< AnalyticEtableEnergy >( src ),
	intrares_evaluator_( src.intrares_evaluator_ ),
	interres_evaluator_( src.interres_evaluator_ )
{
}

methods::EnergyMethodOP
AnalyticEtableEnergy::clone() const {
	return methods::EnergyMethodOP( new AnalyticEtableEnergy( *this ) );
}

void
AnalyticEtableEnergy::setup_for_scoring_( pose::Pose const &pose, scoring::ScoreFunction const& ) const
{
	if ( pose.total_residue() ) {
		if ( pose.residue(1).type().atom_type_set_ptr() != etable().atom_set().lock() ) {
			utility_exit_with_message( "Illegal attempt to score with non-identical atom set between pose and etable " );
		}
	}

	// For debugging if etable options are being updated
	//std::cout << "Check!" << etable().get_lj_hbond_hdis() << " ";
	//std::cout << etable().get_lj_hbond_OH_donor_dis() << std::endl;
}

bool
AnalyticEtableEnergy::defines_intrares_energy(
	EnergyMap const & weights
) const
{
	if ( do_classic_intrares() ) {
		return ( weights[ fa_intra_atr ] != 0 ||
			weights[ fa_intra_rep ] != 0 ||
			weights[ fa_intra_sol ] != 0 );
	}
	if ( put_intra_into_total() ) {
		return ( weights[ fa_atr ] != 0 ||
			weights[ fa_rep ] != 0 ||
			weights[ fa_sol ] != 0 );
	}
	return ( weights[ fa_intra_atr_xover4 ] != 0 ||
		weights[ fa_intra_rep_xover4 ] != 0 ||
		weights[ fa_intra_sol_xover4 ] != 0 );
}

/// @brief
void
AnalyticEtableEnergy::eval_intrares_energy(
	conformation::Residue const & res,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( res.has_variant_type( core::chemical::REPLONLY ) ) {
		return;
	}

	EnergyMap tbemap;
	if ( pose.energies().use_nblist() ) return; // intraresidue atom pairs present in neighborlist, evaluated during finalize

	count_pair::CountPairFunctionOP cpfxn = get_intrares_countpair( res, pose, sfxn );
	//cpfxn->residue_atom_pair_energy( res, res, *this, tbemap );
	intrares_evaluator_.residue_atom_pair_energy( res, res, *cpfxn, tbemap );
	emap[ intrares_evaluator_.st_atr() ] = tbemap[ intrares_evaluator_.st_atr() ];
	emap[ intrares_evaluator_.st_rep() ] = tbemap[ intrares_evaluator_.st_rep() ];
	emap[ intrares_evaluator_.st_sol() ] = tbemap[ intrares_evaluator_.st_sol() ];
}
core::Size
AnalyticEtableEnergy::version() const
{
	return 2; // apl - 2012/06/14 -- analytic version kept in sync with table-based version.
	// return 1; // Initial versioning
}


EtableEvaluator::EtableEvaluator( Etable const & etable ) :
	atr_weight_( 1.0 ),
	rep_weight_( 1.0 ),
	sol_weight_( 1.0 ),
	st_atr_( fa_atr),
	st_rep_( fa_rep ),
	st_sol_( fa_sol ),
	hydrogen_interaction_cutoff2_( etable.hydrogen_interaction_cutoff2() )
{}

EtableEvaluator::~EtableEvaluator() {}

TableLookupEvaluator::TableLookupEvaluator(
	Etable const & etable_in
) :
	EtableEvaluator( etable_in ),
	ljatr_( etable_in.ljatr() ),
	ljrep_( etable_in.ljrep() ),
	solv1_( etable_in.solv1() ),
	solv2_( etable_in.solv2() ),
	dljatr_( etable_in.dljatr() ),
	dljrep_( etable_in.dljrep() ),
	dsolv_( etable_in.dsolv() ),
	safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	etable_bins_per_A2_( etable_in.get_bins_per_A2() ),
	dis2_step_( 1.0 / (Real) etable_bins_per_A2_ )
{}

TableLookupEvaluator::~TableLookupEvaluator() {}


/// @details atom-pair-energy inline type resolution function
void
TableLookupEvaluator::residue_atom_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
TableLookupEvaluator::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_backbone( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
TableLookupEvaluator::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_whole( rsd1, rsd2, *this, emap );
}

//    void
//    TableLookupEvaluator::atom_pair_lk_energy_and_deriv_v_efficient(
//                                                                    conformation::Atom const & atom1,
//                                                                    conformation::Atom const & atom2,
//                                                                    Real & solv1,
//                                                                    Real & dsolv1,
//                                                                    bool const eval_deriv /* = false */
//                                                                    ) const
//    {
//
//        int disbin;
//        Real frac, d2;
//
//        if (interpolate_bins(atom1,atom2,d2,disbin,frac)) {
//
//            int const l1 = solv1_.index( disbin, atom2.type(), atom1.type()),
//            l2 = l1 + 1;
//
//            Real const e1 = solv1_[ l1 ];
//            solv1 = ( e1 + frac * ( solv1_[ l2 ] - e1 ) );
//
//            if ( eval_deriv ){
//                // Following (commented out) is used in dE_dR_over_R below,
//                //  but its a mistake, I think -- rhiju.
//                //   Real e1 = dsolv1_[ l1 ];
//                //   deriv = ( e1 + frac * ( dsolv1_[ l2 ] - e1 ) );
//                dsolv1 = ( solv1_[ l2 ] - solv1_[ l1 ] ) * etable_bins_per_A2_ * std::sqrt( d2 ) * 2;
//            }
//
//        } //if within cutoff
//
//    }


/// @details atom-pair-energy inline type resolution function
void
TableLookupEvaluator::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_backbone_backbone( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
TableLookupEvaluator::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_sidechain( rsd1, rsd2, *this, emap );
}


/// @details first level polymorphic type resolution function
void
TableLookupEvaluator::trie_vs_trie(
	trie::RotamerTrieBase const & trie1,
	trie::RotamerTrieBase const & trie2,
	trie::TrieCountPairBase & cp,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
) const
{
	trie1.trie_vs_trie( trie2, cp, *this, pair_energy_table, temp_table );
}


/// @details first level polymorphic type resolution function
void
TableLookupEvaluator::trie_vs_path(
	trie::RotamerTrieBase const & trie1,
	trie::RotamerTrieBase const & trie2,
	trie::TrieCountPairBase & cp,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
) const
{
	trie1.trie_vs_path( trie2, cp, *this, pair_energy_vector, temp_vector );
}


AnalyticEtableEvaluator::AnalyticEtableEvaluator( Etable const & etable ) :
	EtableEvaluator( etable ),
	etable_( etable ),
	safe_max_dis2_( etable.get_safe_max_dis2() )
{}

AnalyticEtableEvaluator::~AnalyticEtableEvaluator() {}


/// @details atom-pair-energy inline type resolution function
void
AnalyticEtableEvaluator::residue_atom_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
AnalyticEtableEvaluator::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_backbone( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
AnalyticEtableEvaluator::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_whole( rsd1, rsd2, *this, emap );
}


/// @details atom-pair-energy inline type resolution function
void
AnalyticEtableEvaluator::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_backbone_backbone( rsd1, rsd2, *this, emap );
}

/// @details atom-pair-energy inline type resolution function
void
AnalyticEtableEvaluator::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	count_pair::CountPairFunction const & cp,
	EnergyMap & emap
) const
{
	cp.residue_atom_pair_energy_sidechain_sidechain( rsd1, rsd2, *this, emap );
}


/// @details first level polymorphic type resolution function
void
AnalyticEtableEvaluator::trie_vs_trie(
	trie::RotamerTrieBase const & trie1,
	trie::RotamerTrieBase const & trie2,
	trie::TrieCountPairBase & cp,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
) const
{
	trie1.trie_vs_trie( trie2, cp, *this, pair_energy_table, temp_table );
}

/// @details first level polymorphic type resolution function
void
AnalyticEtableEvaluator::trie_vs_path(
	trie::RotamerTrieBase const & trie1,
	trie::RotamerTrieBase const & trie2,
	trie::TrieCountPairBase & cp,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
) const
{
	trie1.trie_vs_path( trie2, cp, *this, pair_energy_vector, temp_vector );
}

void
AnalyticEtableEvaluator::atom_pair_lk_energy_and_deriv_v(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & solE1,
	Real & dsolE1,
	bool const eval_deriv
) const
{
	// could save time by not computing solE2.
	Real solE2;
	etable_.analytic_lk_energy( atom1, atom2, solE1, solE2 );

	if ( eval_deriv ) {
		// could save time by not computing dsolE2.
		Real dsolE2, inv_d;
		etable_.analytic_lk_derivatives( atom1, atom2, dsolE1, dsolE2, inv_d );
	}
}

void
AnalyticEtableEvaluator::atom_pair_lk_energy_and_deriv_v_efficient(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & solE1,
	Real & solE2,
	Real & dsolE1,
	bool const eval_deriv
) const
{
	// could save time by not computing solE2.
	etable_.analytic_lk_energy( atom1, atom2, solE1, solE2 );

	if ( eval_deriv ) {
		// could save time by not computing dsolE2.
		Real dsolE2, inv_d;
		etable_.analytic_lk_derivatives( atom1, atom2, dsolE1, dsolE2, inv_d );
	}
}


} // namespace etable
} // namespace scoring
} // namespace core
