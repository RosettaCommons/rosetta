// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/elec/RNA_FA_ElecEnergy.cc
/// @brief  Electrostatics energy method for RNA class implementation
/// @author Rhiju Das
/// Edited by Joseph Yesselm (9.6.13)


// Unit headers
#include <core/scoring/elec/RNA_FA_ElecEnergy.hh>
#include <core/scoring/elec/RNA_FA_ElecEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/RotamerSetBase.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


// ObjexxFCL headers


// C++

/////////////////////////////////////////////////////////////////////////////////////////
///
/// Hacky (hence the name) implementation of 10r dielectric model, cutoff at 5.5A
///

//     alternatives: WARSHEL (from ligand.cc)
//     E = 322.0637*q1*q2/r/e(r)
//     if ( r < 3 ) e(r) = 16.55
//     else         e(r) = 1 + 60*(1-exp(-0.1*r))
//     Warshel, A. Russell, S. T., Q. Rev. Biophys., 1984, 17, 283-422
//
//
//     sigmoidal dielectric: (hingerty 1985)
//
//     e(r) = D - (D-D0)/2 [ (rS)**2 + 2rS + 2]exp(-rS)
//     with eg:
//     D = 78, D0 = 1, S = 0.3565 (rouzina&bloomfield)
//     D = 80, D0 = 4, S = 0.4 (rohs)

namespace core {
namespace scoring {
namespace elec {


/// @details This must return a fresh instance of the RNA_FA_ElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_FA_ElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new RNA_FA_ElecEnergy( options ) );
}

ScoreTypes
RNA_FA_ElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_elec_rna_phos_phos );
	sts.push_back( fa_elec_rna_phos_sugr );
	sts.push_back( fa_elec_rna_phos_base );
	sts.push_back( fa_elec_rna_sugr_sugr );
	sts.push_back( fa_elec_rna_sugr_base );
	sts.push_back( fa_elec_rna_base_base );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
RNA_FA_ElecEnergy::RNA_FA_ElecEnergy(
	methods::EnergyMethodOptions const & options
):
	parent( options )
{
	set_score_types( methods::EnergyMethodCreatorOP( new RNA_FA_ElecEnergyCreator ) );
}


////////////////////////////////////////////////////////////////////////////
RNA_FA_ElecEnergy::RNA_FA_ElecEnergy( RNA_FA_ElecEnergy const & src ):
	parent( src )
{
	set_score_types( methods::EnergyMethodCreatorOP( new RNA_FA_ElecEnergyCreator ) );
}

/// clone
methods::EnergyMethodOP
RNA_FA_ElecEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_FA_ElecEnergy( *this ) );
}


void
RNA_FA_ElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const &) const
{
	pose.update_residue_neighbors();

}


void
RNA_FA_ElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


// The FA_ElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
RNA_FA_ElecEnergy::setup_for_packing(
	pose::Pose &,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
RNA_FA_ElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
RNA_FA_ElecEnergy::update_residue_for_packing(
	pose::Pose &,
	Size
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

RNAAtomType
assign_rna_atom_type(
	conformation::Residue const & rsd,
	Size const i
) {
	if ( rsd.atom_is_hydrogen(i) ) {
		return assign_rna_atom_type( rsd, rsd.atom_base( i ) );
	}

	// AMW TODO: this cannot stand.
	if ( i == 1 || i == 2 || i == 3 || i == 4 || i == 9 ) { //P, OP2, OP1, O5', O3'
		return PHOSPHATE;
	} else if ( i < rsd.first_sidechain_atom() || i == 12 ) {
		return SUGAR;
	} else {
		return BASE;
	}
}


///
void
RNA_FA_ElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const & scfxn,
	EnergyMap & emap
) const {
	if ( ! rsd1.is_RNA() || ! rsd2.is_RNA() ) { return; }

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_phos) ) {
		emap[ fa_elec_rna_phos_phos ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_sugr_sugr) ) {
		emap[ fa_elec_rna_sugr_sugr ] += rna_fa_elec_one_way(rsd1,rsd2,SUGAR,SUGAR);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_base_base) ) {
		emap[ fa_elec_rna_base_base ] += rna_fa_elec_one_way(rsd1,rsd2,BASE,BASE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_sugr) ) {
		emap[ fa_elec_rna_phos_sugr ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,SUGAR) +
			rna_fa_elec_one_way(rsd1,rsd2,SUGAR,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_phos_base) ) {
		emap[ fa_elec_rna_phos_base ] += rna_fa_elec_one_way(rsd1,rsd2,PHOSPHATE,BASE) +
			rna_fa_elec_one_way(rsd1,rsd2,BASE,PHOSPHATE);
	}

	if ( scfxn.has_nonzero_weight(fa_elec_rna_sugr_base) ) {
		emap[ fa_elec_rna_sugr_base ] += rna_fa_elec_one_way(rsd1,rsd2,SUGAR,BASE) +
			rna_fa_elec_one_way(rsd1,rsd2,BASE,SUGAR);
	}
}


//////////////////////////////////////////////////////////////////////////////////
// Sept 5th 2013, implemetnation of an optmized design function only phos-phos interactions are evaluated!
//////////////////////////////////////////////////////////////////////////////////




void
RNA_FA_ElecEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	if ( set1.num_rotamers() >= 1 && set2.num_rotamers() >= 1 &&
			set1.rotamer(1)->is_RNA() && set2.rotamer(1)->is_RNA() ) {
		grandparent::evaluate_rotamer_pair_energies( set1, set2, pose, sfxn, weights, energy_table );
	} // else, non rna rna interaction; early return
}

void
RNA_FA_ElecEnergy::backbone_backbone_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	emap[ fa_elec_rna_phos_phos ] +=  rna_fa_elec_one_way(rsd1, rsd2, PHOSPHATE, PHOSPHATE);
}


Real
RNA_FA_ElecEnergy::rna_fa_elec_one_way(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	RNAAtomType const & type1 ,
	RNAAtomType const & type2
) const {
	Real energy( 0.0 );

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size i=1, i_end = rsd1.natoms(); i<= i_end; ++i ) {

		RNAAtomType const atom_type_1 = assign_rna_atom_type(rsd1, i);
		if ( atom_type_1 != type1 ) continue;

		Vector const & i_xyz( rsd1.xyz(i) );
		Real const i_charge( rsd1.atomic_charge(i) );

		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {

			RNAAtomType const atom_type_2 = assign_rna_atom_type(rsd2, j);
			if ( atom_type_2 != type2 ) continue;

			Real weight(1.0);
			Size path_dist( 0 );
			if ( ! cpfxn->count( i, j, weight, path_dist ) ) continue;

			energy += weight * coulomb().eval_atom_atom_fa_elecE( i_xyz, i_charge, rsd2.xyz(j), rsd2.atomic_charge(j));
		}
	}
	return energy;
}

void
RNA_FA_ElecEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	if ( set.num_rotamers() >= 1 && set.rotamer(1)->is_RNA() && residue.is_RNA() ) {
		grandparent::evaluate_rotamer_background_energies( set, residue, pose, sfxn, weights, energy_vector );
	} // else, non rna rna interaction; early return
}



//////////////////////////////////////
bool is_phosphate_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_phosphate_2( rsd, rsd.atom_base( i ) );
	} else {
		//MAGIC NUMBERS! BAD!
		return ( i==1 || i==2 || i==3 || i==4 || i==9 ); //P, OP2, OP1, O5', O3'
	}
}

//////////////////////////////////////
bool is_sugar_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_sugar_2( rsd, rsd.atom_base( i ) );
	} else {
		return ( ( i < rsd.first_sidechain_atom() && !is_phosphate_2(rsd, i) ) ||
			i == 12 /*hey this is really bad!*/ );
	}
}

//////////////////////////////////////
bool is_base_2( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_base_2( rsd, rsd.atom_base( i ) );
	} else {
		return (!is_sugar_2(rsd, i) && !is_phosphate_2(rsd, i) );
	}
}

void
RNA_FA_ElecEnergy::eval_atom_derivative_RNA(
	conformation::Residue const & rsd1,
	Size const & i,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real const i_charge( rsd1.atomic_charge( i ) );
	Vector const & i_xyz( rsd1.xyz(i) );

	// NOTE THAT is_base, is_sugar, and is_phosphate contains MAGIC NUMBERS (for speed!)
	// Might be better to make it precomputed as part of the residue_type definition?
	//This repeats some unnecessary stuff, but hey, derivatives don't have to be that fast.
	bool const atom1_is_base = is_base_2(rsd1, i);
	bool const atom1_is_sugar = is_sugar_2(rsd1, i);
	bool const atom1_is_phosphate = is_phosphate_2(rsd1, i);
	debug_assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

	for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {
		Real const j_charge( rsd2.atomic_charge(j) );
		if ( j_charge == 0.0 ) continue;
		Real weight(1.0);
		Size path_dist( 0 );
		if ( ! cpfxn->count( i, j, weight, path_dist ) ) continue;

		Vector const & j_xyz( rsd2.xyz(j) );
		Vector const f2( i_xyz - j_xyz );
		Real const dis2( f2.length_squared() );
		Real const dE_dr_over_r = weight *
			coulomb().eval_dfa_elecE_dr_over_r( dis2, i_charge, j_charge );
		if ( dE_dr_over_r == 0.0 ) continue;

		Vector const f1( i_xyz.cross( j_xyz ) );

		bool const atom2_is_base = is_base_2(rsd2, j);
		bool const atom2_is_sugar = is_sugar_2(rsd2, j);
		bool const atom2_is_phosphate = is_phosphate_2(rsd2, j);
		debug_assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

		if ( atom1_is_base && atom2_is_base ) {
			F1 += weights[ fa_elec_rna_base_base ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_base_base ] * dE_dr_over_r * f2;
		} else if (  (atom1_is_base && atom2_is_sugar)  || (atom1_is_sugar && atom2_is_base) ) {
			F1 += weights[ fa_elec_rna_sugr_base ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_sugr_base ] * dE_dr_over_r * f2;
		} else if (  (atom1_is_base && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_base) ) {
			F1 += weights[ fa_elec_rna_phos_base ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_phos_base ] * dE_dr_over_r * f2;
		} else if (  (atom1_is_sugar && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_sugar) ) {
			F1 += weights[ fa_elec_rna_phos_sugr ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_phos_sugr ] * dE_dr_over_r * f2;
		} else if (  (atom1_is_sugar && atom2_is_sugar)  ) {
			F1 += weights[ fa_elec_rna_sugr_sugr ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_sugr_sugr ] * dE_dr_over_r * f2;
		} else if (  (atom1_is_phosphate && atom2_is_phosphate)  ) {
			F1 += weights[ fa_elec_rna_phos_phos ] * dE_dr_over_r * f1;
			F2 += weights[ fa_elec_rna_phos_phos ] * dE_dr_over_r * f2;
		}
	}
}


void
RNA_FA_ElecEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace etable::count_pair;

	// what is my charge?
	Size const pos1( atom_id.rsd() );
	Size const i   ( atom_id.atomno() );
	conformation::Residue const & rsd1( pose.residue( pos1 ) );

	if ( ! rsd1.is_RNA() ) return;

	Real const i_charge( rsd1.atomic_charge( i ) );
	int const pos1_map( domain_map( pos1 ) );
	bool const pos1_fixed( pos1_map != 0 );

	if ( i_charge == 0.0 ) return;

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const pos2( (*iru)->get_other_ind( pos1 ) );
		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );
		debug_assert( pos2 != pos1 );

		if ( rsd2.is_RNA() ) {
			eval_atom_derivative_RNA( rsd1, i, rsd2, weights, F1, F2 );
		}
	} // loop over nbrs of rsd1
}


/// @brief RNA_FA_ElecEnergy2 is context independent; no context graphs required
void
RNA_FA_ElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
RNA_FA_ElecEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
