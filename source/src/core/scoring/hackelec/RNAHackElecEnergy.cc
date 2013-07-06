// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hackelec/RNAHackElecEnergy.cc
/// @brief  Electrostatics energy method for RNA class implementation
/// @author Rhiju Das


// Unit headers
#include <core/scoring/hackelec/RNAHackElecEnergy.hh>
#include <core/scoring/hackelec/RNAHackElecEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/RotamerSetBase.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSetFactory.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


// ObjexxFCL headers



// C++

/////////////////////////////////////////////////////////////////////////////////////////
///
/// Hacky (hence the name) implementation of 10r dielectric model, cutoff at 5.5A
///
///

//
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
namespace hackelec {


/// @details This must return a fresh instance of the RNAHackElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNAHackElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new RNAHackElecEnergy( options );
}

ScoreTypes
RNAHackElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hack_elec_rna_phos_phos );
	sts.push_back( hack_elec_rna_phos_sugr );
	sts.push_back( hack_elec_rna_phos_base );
	sts.push_back( hack_elec_rna_sugr_sugr );
	sts.push_back( hack_elec_rna_sugr_base );
	sts.push_back( hack_elec_rna_base_base );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
RNAHackElecEnergy::RNAHackElecEnergy(
	methods::EnergyMethodOptions const & options
):
	parent( options )
{
	set_score_types( new RNAHackElecEnergyCreator );
}


////////////////////////////////////////////////////////////////////////////
RNAHackElecEnergy::RNAHackElecEnergy( RNAHackElecEnergy const & src ):
	parent( src )
{
	set_score_types( new RNAHackElecEnergyCreator );
}

/// clone
methods::EnergyMethodOP
RNAHackElecEnergy::clone() const
{
	return new RNAHackElecEnergy( *this );
}

///
void
RNAHackElecEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

///
void
RNAHackElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}


// The HackElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
RNAHackElecEnergy::setup_for_packing(
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
RNAHackElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase &
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}


// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
RNAHackElecEnergy::update_residue_for_packing(
	pose::Pose &,
	Size
) const
{
	/// noop -- this noop is essential for overriding the base class behavior.
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
RNAHackElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( rsd1.is_RNA() && rsd2.is_RNA() ) {
		residue_pair_energy_RNA( rsd1, rsd2, emap ); // In the grand spirit of hack_elec hackiness
	}
	//std::cout << rsd1.seqpos() << ' ' << rsd2.seqpos() << ' ' << score << std::endl;
}


//////////////////////////////////////
bool is_phosphate( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_phosphate( rsd, rsd.atom_base( i ) );
	} else {
		//MAGIC NUMBERS! BAD!
		return ( i==1 || i==2 || i==3 || i==4 || i==9 ); //P, O1P, O2P, O5*, O3*
	}
}

//////////////////////////////////////
bool is_sugar( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_sugar( rsd, rsd.atom_base( i ) );
	} else {
		return ( ( i < rsd.first_sidechain_atom() && !is_phosphate(rsd, i) ) ||
						 i == 12 /*hey this is really bad!*/ );
	}
}

//////////////////////////////////////
bool is_base( conformation::Residue const & rsd, Size const i )
{
	if ( rsd.atom_is_hydrogen( i ) ) {
		return is_base( rsd, rsd.atom_base( i ) );
	} else {
		return (!is_sugar(rsd, i) && !is_phosphate(rsd, i) );
	}
}

//////////////////////////////////////////////////////////////////////////////////
// July 22 2011. OK could potentially optimize by calculating the energy for the atom-atom pairs in which the weight term is non-zero.
// For example, for standard weight only hack_elec_rna_phos_phos weight is non-zero. So for this case can skip non-phosphate atoms.
// Following copies a little code, but at least separates out this inelegant and
// probably useless RNA stuff from the usual stuff.
Real
RNAHackElecEnergy::residue_pair_energy_RNA(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	EnergyMap & emap
) const
{
	assert( rsd1.is_RNA() );
	assert( rsd2.is_RNA() );

	using namespace etable::count_pair;

	Real total_score( 0.0 );

	CountPairFunctionOP cpfxn =
		CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	for ( Size i=1, i_end = rsd1.natoms(); i<= i_end; ++i ) {
		Vector const & i_xyz( rsd1.xyz(i) );
		Real const i_charge( rsd1.atomic_charge(i) );

		// NOTE THAT is_base, is_sugar, and is_phosphate contains MAGIC NUMBERS (for speed!)
		// Might be better to make it precomputed as part of the residue_type definition?
		bool const atom1_is_base = is_base(rsd1, i);
		bool const atom1_is_sugar = is_sugar(rsd1, i);
		bool const atom1_is_phosphate = is_phosphate(rsd1, i);
		assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

		//		if (atom1_is_base) std::cout << "THIS BETTER BE A BASE ATOM: "  << rsd1.atom_name( i ) << std::endl;

		if ( i_charge == 0.0 ) continue;
		for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {
			Real const j_charge( rsd2.atomic_charge(j) );
			if ( j_charge == 0.0 ) continue;
			Real weight(1.0);
			Size path_dist( 0 );
			if ( cpfxn->count( i, j, weight, path_dist ) ) {
				Real score = weight *
					coulomb().eval_atom_atom_hack_elecE( i_xyz, i_charge, rsd2.xyz(j), j_charge);

				total_score += score;

				bool const atom2_is_base = is_base(rsd2, j);
				bool const atom2_is_sugar = is_sugar(rsd2, j);
				bool const atom2_is_phosphate = is_phosphate(rsd2, j);
				assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

				if ( atom1_is_base && atom2_is_base ) {
					emap[ hack_elec_rna_base_base ] += score;
				} else if (  (atom1_is_base && atom2_is_sugar)  || (atom1_is_sugar && atom2_is_base) ) {
					emap[ hack_elec_rna_sugr_base ] += score;
				} else if (  (atom1_is_base && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_base) ) {
					emap[ hack_elec_rna_phos_base ] += score;
				} else if (  (atom1_is_sugar && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_sugar) ) {
					emap[ hack_elec_rna_phos_sugr ] += score;
				} else if (  (atom1_is_sugar && atom2_is_sugar)  ) {
					emap[ hack_elec_rna_sugr_sugr ] += score;
				} else if (  (atom1_is_phosphate && atom2_is_phosphate)  ) {
					emap[ hack_elec_rna_phos_phos ] += score;
				} else {
					std::cout << "PROBLEM! " << rsd1.atom_name( i ) << " " << rsd2.atom_name( j ) << std::endl;
				}
			}
		}
	}

	return total_score;

}

void
RNAHackElecEnergy::evaluate_rotamer_pair_energies(
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
RNAHackElecEnergy::evaluate_rotamer_background_energies(
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


void
RNAHackElecEnergy::eval_atom_derivative(
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

	//Vector const & i_xyz( rsd1.xyz(i) );

	// cached energies object
	Energies const & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

// 	kinematics::DomainMap const & domain_map( energies.domain_map() );
// 	bool const pos1_fixed( !energies.res_moved( pos1 ) );
// 	assert( pos1_fixed == ( domain_map(pos1) != 0 ) ); // this is probably not generally true but I'm curious

	// loop over *all* nbrs of rsd1 (not just upper or lower)
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( pos1 )->const_edge_list_begin(),
			irue = energy_graph.get_node( pos1 )->const_edge_list_end();
			iru != irue; ++iru ) {
		//EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		//Size const pos2( edge->get_second_node_ind() );
		Size const pos2( (*iru)->get_other_ind( pos1 ) );

		if ( pos1_fixed && pos1_map == domain_map( pos2 ) ) continue; // fixed wrt one another

		conformation::Residue const & rsd2( pose.residue( pos2 ) );

		assert( pos2 != pos1 );

		if ( rsd2.is_RNA() ) {
			eval_atom_derivative_RNA( rsd1, i, rsd2, weights, F1, F2 );
		}

	} // loop over nbrs of rsd1

}


//////////////////////////////////////////////////
void
RNAHackElecEnergy::eval_atom_derivative_RNA(
	conformation::Residue const & rsd1,
	Size const & i,
	conformation::Residue const & rsd2,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
 	) const
{

	using namespace etable::count_pair;

	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	Real const i_charge( rsd1.atomic_charge( i ) );
	Vector const & i_xyz( rsd1.xyz(i) );


	// NOTE THAT is_base, is_sugar, and is_phosphate contains MAGIC NUMBERS (for speed!)
	// Might be better to make it precomputed as part of the residue_type definition?
	//This repeats some unnecessary stuff, but hey, derivatives don't have to be that fast.
	bool const atom1_is_base = is_base(rsd1, i);
	bool const atom1_is_sugar = is_sugar(rsd1, i);
	bool const atom1_is_phosphate = is_phosphate(rsd1, i);
	assert( atom1_is_base || atom1_is_sugar || atom1_is_phosphate );

	for ( Size j=1, j_end = rsd2.natoms(); j<= j_end; ++j ) {
		Real const j_charge( rsd2.atomic_charge(j) );
		if ( j_charge == 0.0 ) continue;
		Real weight(1.0);
		Size path_dist( 0 );
		if ( cpfxn->count( i, j, weight, path_dist ) ) {
			Vector const & j_xyz( rsd2.xyz(j) );
			Vector const f2( i_xyz - j_xyz );
			Real const dis2( f2.length_squared() );
			Real const dE_dr_over_r = weight *
				coulomb().eval_dhack_elecE_dr_over_r( dis2, i_charge, j_charge );
			if ( dE_dr_over_r != 0.0 ) {
				Vector const f1( i_xyz.cross( j_xyz ) );
				//F1 += weights[ hack_elec ] * dE_dr_over_r * f1;
				//F2 += weights[ hack_elec ] * dE_dr_over_r * f2;

				bool const atom2_is_base = is_base(rsd2, j);
				bool const atom2_is_sugar = is_sugar(rsd2, j);
				bool const atom2_is_phosphate = is_phosphate(rsd2, j);
				assert( atom2_is_base || atom2_is_sugar || atom2_is_phosphate );

				if ( atom1_is_base && atom2_is_base ) {
					F1 += weights[ hack_elec_rna_base_base ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_base_base ] * dE_dr_over_r * f2;
				} else if (  (atom1_is_base && atom2_is_sugar)  || (atom1_is_sugar && atom2_is_base) ) {
					F1 += weights[ hack_elec_rna_sugr_base ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_sugr_base ] * dE_dr_over_r * f2;
				} else if (  (atom1_is_base && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_base) ) {
					F1 += weights[ hack_elec_rna_phos_base ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_phos_base ] * dE_dr_over_r * f2;
				} else if (  (atom1_is_sugar && atom2_is_phosphate)  || (atom1_is_phosphate && atom2_is_sugar) ) {
					F1 += weights[ hack_elec_rna_phos_sugr ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_phos_sugr ] * dE_dr_over_r * f2;
				} else if (  (atom1_is_sugar && atom2_is_sugar)  ) {
					F1 += weights[ hack_elec_rna_sugr_sugr ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_sugr_sugr ] * dE_dr_over_r * f2;
				} else if (  (atom1_is_phosphate && atom2_is_phosphate)  ) {
					F1 += weights[ hack_elec_rna_phos_phos ] * dE_dr_over_r * f1;
					F2 += weights[ hack_elec_rna_phos_phos ] * dE_dr_over_r * f2;
				}

			}
		}
	}
}


/// @brief RNAHackElecEnergy is context independent; no context graphs required
void
RNAHackElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
RNAHackElecEnergy::version() const
{
	return 1; // Initial versioning
}



}
}
}
