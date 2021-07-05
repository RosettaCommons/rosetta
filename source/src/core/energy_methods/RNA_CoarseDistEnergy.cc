// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/energy_methods/RNA_CoarseDistEnergy.cc
/// @brief Score two-body energies in coarse RNA poses between P, S, and CEN using a statistical potential.
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <core/energy_methods/RNA_CoarseDistEnergy.hh>

#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/OneDDistPotential.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/MinimizationData.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/NeighborList.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

#include <core/energy_methods/RNA_CoarseDistEnergyCreator.hh> // AUTO IWYU For RNA_CoarseDistEnergyCreator

static basic::Tracer TR( "core.energy_methods.RNA_CoarseDistEnergy" );

namespace core {
namespace energy_methods {

using namespace core::scoring;

RNA_CoarseDistEnergy::RNA_CoarseDistEnergy( core::scoring::methods::EnergyMethodOptions const & /*options*/ ):
	core::scoring::methods::ContextIndependentTwoBodyEnergy( utility::pointer::make_shared< RNA_CoarseDistEnergyCreator >() ),
	P_P_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_P_P.bin" ) ) ),
	P_S_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_P_S.bin" ) ) ),
	P_CEN_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_P_CEN.bin" ) ) ),
	S_S_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_S_S.bin" ) ) ),
	S_CEN_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_S_CEN.bin" ) ) ),
	CEN_CEN_potential_( core::scoring::ScoringManager::get_instance()->get_OneDDistPotential( basic::database::full_name("scoring/rna/coarse/dist_CEN_CEN.bin" ) ) )
{
}

RNA_CoarseDistEnergy::~RNA_CoarseDistEnergy(){}

/// clone
methods::EnergyMethodOP
RNA_CoarseDistEnergy::clone() const
{
	return utility::pointer::make_shared< RNA_CoarseDistEnergy >( *this );
}


inline
Real
RNA_CoarseDistEnergy::score_atom_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Size const at1,
	Size const at2
) const
{
	Real energy = 0;
	if ( at1 == 1 && at2 == 1 ) { // sum 2
		energy = P_P_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining P-P potential for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  P_P_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	} else if ( (at1 == 1 && at2 == 2) || (at1 == 2 && at2 == 1) ) { // sum 3
		energy = P_S_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining P_S_potential_ for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  P_S_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	} else if ( (at1 == 1 && at2 == 3) || (at1 == 3 && at2 == 1) ) { // sum 4
		energy = P_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining P_CEN_potential_ for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  P_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	} else if ( (at1 == 2 && at2 == 2) ) { // aw, sum 4
		energy = S_S_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining S_S_potential_ for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  S_S_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	} else if ( ( at1 == 2 && at2 == 3) || (at1 == 3 && at2 == 2) ) {
		energy = S_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining S_CEN_potential_ for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  S_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	} else if ( at1 == 3 && at2 == 3 ) {
		energy = CEN_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
		if ( energy < -6 ) {
			TR << "Examining CEN_CEN_potential_ for " << rsd1.seqpos() << " " << rsd2.seqpos() << ". Distance " << rsd1.xyz(at1).distance(rsd2.xyz(at2)) << std::endl;
			TR << "Value of " <<  CEN_CEN_potential_.evaluate( rsd1.xyz(at1).distance(rsd2.xyz(at2))) << std::endl;
		}
	}

	return energy;
}

inline
Real
RNA_CoarseDistEnergy::deriv_atom_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Size const at1,
	Size const at2
) const
{
	Real energy = 0;
	if ( at1 == 1 && at2 == 1 ) { // sum 2
		energy = P_P_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	} else if ( (at1 == 1 && at2 == 2) || (at1 == 2 && at2 == 1) ) { // sum 3
		energy = P_S_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	} else if ( (at1 == 1 && at2 == 3) || (at1 == 3 && at2 == 1) ) { // sum 4
		energy = P_CEN_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	} else if ( at1 == 2 && at2 == 2 ) { // aw, sum 4
		energy = S_S_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	} else if ( (at1 == 2 && at2 == 3) || (at1 == 3 && at2 == 2 ) ) {
		energy = S_CEN_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	} else if ( at1 == 3 && at2 == 3 ) {
		energy = CEN_CEN_potential_.get_derivative( rsd1.xyz(at1).distance(rsd2.xyz(at2)));
	}

	return energy;
}

void
RNA_CoarseDistEnergy::residue_pair_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & ,//sfxn,
	EnergyMap & emap
) const {
	// AMW TODO: assert we're coarse_rna
	if ( pose.energies().use_nblist() ) return;

	Real score(0.0);

	if ( ! defines_score_for_residue_pair( rsd1, rsd2, true ) ) return;
	// TR << "good on " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

	// NULL if no info
	if ( !rsd1.is_bonded( rsd2 ) && !rsd1.is_pseudo_bonded( rsd2 ) ) {
		for ( Size ii = 1; ii <= 3; ++ii ) {
			for ( Size jj = 1; jj <= 3; ++jj ) {
				// TR << "About to evaluate for " << ii << jj << std::endl;
				score += score_atom_pair( rsd1, rsd2, ii, jj );
			}
		}
	}
	emap[ rna_coarse_dist ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

bool
RNA_CoarseDistEnergy::defines_score_for_residue_pair(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const {

	if ( res1.seqpos() == res2.seqpos() ) {
		return false;
	}

	return res_moving_wrt_eachother;
}

bool
RNA_CoarseDistEnergy::use_extended_residue_pair_energy_interface() const { return true; }

void
RNA_CoarseDistEnergy::residue_pair_energy_ext(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	core::pose::Pose const & pose,
	ScoreFunction const & ,//sfxn,
	EnergyMap & emap
) const {
	if ( pose.energies().use_nblist_auto_update() ) return;

	// reuse the elec pair nblist. it's fine!

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( elec_pair_nblist ) ));
	auto const & nblist( static_cast< core::scoring::ResiduePairNeighborList const & > ( min_data.get_data_ref( elec_pair_nblist ) ) );
	Real score( 0.0 );
	utility::vector1< core::scoring::SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	for ( auto const & neighb : neighbs ) {
		if ( neighb.atomno1() > 3 || neighb.atomno2() > 3 ) continue;
		score += score_atom_pair( rsd1, rsd2, neighb.atomno1(), neighb.atomno2() );
	}
	emap[ rna_coarse_dist ] += score;
}

void
RNA_CoarseDistEnergy::setup_for_minimizing_for_residue(
	core::conformation::Residue const & ,//rsd,
	core::pose::Pose const & ,//pose,
	ScoreFunction const & ,//sfxn,
	core::kinematics::MinimizerMapBase const & ,//minmap,
	basic::datacache::BasicDataCache & ,//residue_data_cache
	ResSingleMinimizationData & //res_data_cache
) const { }


/// @brief Returns a regular count-pair function as opposed to a CountPairRepresentative function
etable::count_pair::CountPairFunctionCOP
RNA_CoarseDistEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace core::scoring::etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return utility::pointer::make_shared< CountPairNone >();

	// if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
	//  return /*count_pair_full_ ?
	//   CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3FULL ) :*/
	//   CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	//  //return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	// }
	return utility::pointer::make_shared< CountPairAll >();

}

void
RNA_CoarseDistEnergy::setup_for_minimizing_for_residue_pair(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & pose,
	ScoreFunction const & ,//sfxn,
	core::kinematics::MinimizerMapBase const & ,//minmap,
	ResSingleMinimizationData const & ,//res1_data_cache,
	ResSingleMinimizationData const & ,//res2_data_cache,
	ResPairMinimizationData & data_cache
) const {

	if ( pose.energies().use_nblist_auto_update() ) return;

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );
	debug_assert( rsd1.seqpos() < rsd2.seqpos() );

	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist(
		utility::pointer::static_pointer_cast< core::scoring::ResiduePairNeighborList > ( data_cache.get_data( elec_pair_nblist ) ));
	if ( ! nblist ) nblist = utility::pointer::make_shared< ResiduePairNeighborList >();

	// Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	static Real const XX2 = std::pow( 16 + 2*0.75, 2 );

	// if ( !use_cp_rep_ ) {
	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );
	// }
	//  else {
	//  // ensure that the resmaps are initialized for this pair
	//  get_countpair_representative_atom(rsd1.type(),1);
	//  get_countpair_representative_atom(rsd2.type(),1);

	//  nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair, cp_rep_map_->get_map(rsd1.type()), cp_rep_map_->get_map(rsd2.type()) );
	// }

	data_cache.set_data( elec_pair_nblist, nblist );
}

bool
RNA_CoarseDistEnergy::requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( core::pose::Pose const & /*pose*/ ) const {
	return true;
}

void
RNA_CoarseDistEnergy::setup_for_scoring_for_residue(
	core::conformation::Residue const & ,//rsd,
	core::pose::Pose const & ,//pose,
	ScoreFunction const & ,//sfxn,
	ResSingleMinimizationData & //min_data
) const { }

bool
RNA_CoarseDistEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( core::pose::Pose const & /*pose*/ ) const {
	return true;
}

void
RNA_CoarseDistEnergy::setup_for_derivatives_for_residue(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data,
	basic::datacache::BasicDataCache & //res_data_cache
) const
{
	setup_for_scoring_for_residue( rsd, pose, sfxn, min_data );
}


// bool
// RNA_CoarseDistEnergy::requires_a_setup_for_scoring_for_residue_pair_opportunity( core::pose::Pose const & pose ) const;

// void
// RNA_CoarseDistEnergy::setup_for_scoring_for_residue_pair(
//  core::conformation::Residue const & rsd1,
//  core::conformation::Residue const & rsd2,
//  ResSingleMinimizationData const & minsingle_data1,
//  ResSingleMinimizationData const & minsingle_data2,
//  core::pose::Pose const & pose,
//  ScoreFunction const & sfxn,
//  ResPairMinimizationData & data_cache
// ) const;

// bool
// RNA_CoarseDistEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( core::pose::Pose const & pose ) const;

// void
// RNA_CoarseDistEnergy::setup_for_derivatives_for_residue_pair(
//  conformation::Residue const & rsd1,
//  conformation::Residue const & rsd2,
//  ResSingleMinimizationData const & minsingle_data1,
//  ResSingleMinimizationData const & minsingle_data2,
//  pose::Pose const & pose,
//  ScoreFunction const & sfxn,
//  ResPairMinimizationData & data_cache
// ) const;

void
RNA_CoarseDistEnergy::eval_residue_pair_derivatives(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	core::pose::Pose const & pose, // provides context
	EnergyMap const & ,//weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	if ( pose.energies().use_nblist_auto_update() ) return;

	debug_assert( rsd1.seqpos() < rsd2.seqpos() );
	debug_assert( utility::pointer::dynamic_pointer_cast< ResiduePairNeighborList const > (min_data.get_data( elec_pair_nblist ) ));

	auto const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( elec_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );

	// weight_triple wtrip;
	// setup_weight_triple( weights, wtrip );
	for ( auto const & neighb : neighbs ) {
		Size at1 = neighb.atomno1();
		Size at2 = neighb.atomno2();
		if ( at1 > 3 || at2 > 3 ) continue;

		Vector const & atom1xyz = rsd1.xyz( at1 );
		Vector const & atom2xyz = rsd2.xyz( at2 );
		Vector f1( 0.0 ), f2 = ( atom1xyz - atom2xyz );

		Real dE_dr_over_r = deriv_atom_pair( rsd1, rsd2, at1, at2 );

		if ( dE_dr_over_r != 0.0 ) {
			f1 = atom1xyz.cross( atom2xyz );
			Vector f1s = dE_dr_over_r * f1;
			Vector f2s = dE_dr_over_r * f2;
			r1_atom_derivs[ at1 ].f1() += f1s;
			r1_atom_derivs[ at1 ].f2() += f2s;
			r2_atom_derivs[ at2 ].f1() -= f1s;
			r2_atom_derivs[ at2 ].f2() -= f2s;
		}
	}
}

void
RNA_CoarseDistEnergy::backbone_backbone_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & ,//pose,
	ScoreFunction const & ,//sfxn,
	EnergyMap & emap
) const {
	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( !rsd1.is_bonded( rsd2 ) && !rsd1.is_pseudo_bonded( rsd2 ) ) {
		// no countpair!
		score += P_P_potential_.evaluate( rsd1.xyz( 1 ).distance( rsd2.xyz( 1 ) ) );
		score += P_S_potential_.evaluate( rsd1.xyz( 1 ).distance( rsd2.xyz( 2 ) ) );
		score += P_S_potential_.evaluate( rsd1.xyz( 2 ).distance( rsd2.xyz( 1 ) ) );
		score += S_S_potential_.evaluate( rsd1.xyz( 2 ).distance( rsd2.xyz( 2 ) ) );
	}
	emap[ rna_coarse_dist ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

void
RNA_CoarseDistEnergy::backbone_sidechain_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & ,//pose,
	ScoreFunction const & ,//sfxn,
	EnergyMap & emap
) const {
	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( !rsd1.is_bonded( rsd2 ) && !rsd1.is_pseudo_bonded( rsd2 ) ) {
		// no countpair!
		score += P_CEN_potential_.evaluate( rsd1.xyz( 1 ).distance( rsd2.xyz( 3 ) ) );
		score += S_CEN_potential_.evaluate( rsd1.xyz( 2 ).distance( rsd2.xyz( 3 ) ) );
		score += P_CEN_potential_.evaluate( rsd1.xyz( 3 ).distance( rsd2.xyz( 1 ) ) );
		score += S_CEN_potential_.evaluate( rsd1.xyz( 3 ).distance( rsd2.xyz( 2 ) ) );
	}
	emap[ rna_coarse_dist ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

void
RNA_CoarseDistEnergy::sidechain_sidechain_energy(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::pose::Pose const & ,//pose,
	ScoreFunction const & ,//sfxn,
	EnergyMap & emap
) const {
	Real score(0.0);

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	if ( !rsd1.is_bonded( rsd2 ) && !rsd1.is_pseudo_bonded( rsd2 ) ) {
		// no countpair!
		score += CEN_CEN_potential_.evaluate( rsd1.xyz( 3 ).distance( rsd2.xyz( 3 ) ) );
	}
	emap[ rna_coarse_dist ] += score;
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}

void
RNA_CoarseDistEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{
	if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) return;

	// add in contributions from the nblist atom-pairs
	core::scoring::NeighborList const & nblist
		( pose.energies().nblist( core::scoring::EnergiesCacheableDataType::ELEC_NBLIST ) );

	nblist.check_domain_map( pose.energies().domain_map() );

	utility::vector1< conformation::Residue const * > resvect;
	resvect.reserve( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		resvect.push_back( & pose.residue( ii ) );
	}

	Real total_score( 0.0 );
	for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
		conformation::Residue const & ires( *resvect[i] );
		for ( Size ii=1, ii_end=ires.natoms(); ii<= ii_end; ++ii ) {
			core::scoring::AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );
			for ( auto const & nbr : nbrs ) {
				Size const  j( nbr.rsd() );
				Size const jj( nbr.atomno() );

				conformation::Residue const & jres( *resvect[j] );

				Real score = score_atom_pair( ires, jres, ii, jj );
				total_score += score;
			}
		}
	}
	totals[ rna_coarse_dist ] += total_score;
}

// void
// RNA_CoarseDistEnergy::evaluate_rotamer_pair_energies(
//  core::conformation::RotamerSetBase const & set1,
//  core::conformation::RotamerSetBase const & set2,
//  core::pose::Pose const & pose,
//  ScoreFunction const & sfxn,
//  EnergyMap const & weights,
//  ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
// ) const

// void
// RNA_CoarseDistEnergy::evaluate_rotamer_background_energies(
//  core::conformation::RotamerSetBase const & set,
//  core::conformation::Residue const & residue,
//  core::pose::Pose const & pose,
//  ScoreFunction const & sfxn,
//  EnergyMap const & weights,
//  utility::vector1< core::PackerEnergy > & energy_vector
// ) const;

// void
// RNA_CoarseDistEnergy::evaluate_rotamer_background_energy_maps(
//  core::conformation::RotamerSetBase const & set,
//  core::conformation::Residue const & residue,
//  core::pose::Pose const & pose,
//  ScoreFunction const & sfxn,
//  EnergyMap const & weights,
//  utility::vector1< EnergyMap > & emaps
// ) const;


/// @brief Instantiate a new RNA_CoarseDistEnergy
methods::EnergyMethodOP
RNA_CoarseDistEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< RNA_CoarseDistEnergy >( options );
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
ScoreTypes
RNA_CoarseDistEnergyCreator::score_types_for_method() const {
	return { core::scoring::rna_coarse_dist };
}

} //energy_methods
} //core

