// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairIntraResC4.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 4 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#include <core/scoring/etable/count_pair/CountPairGeneric.hh>

#include <core/types.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/PseudoBond.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>


namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

CountPairGeneric::CountPairGeneric(
	conformation::Residue const & res1,
	conformation::Residue const & res2
) :
	//r1_( res1 ),
	//r2_( res2 ),
	n_connect_( 0 ),
	n_pconnect_( 0 ),
	crossover_( 3 )
{
	using namespace conformation;

	debug_assert( res1.seqpos() != res2.seqpos() );

	if ( res1.is_bonded( res2 ) ) {
		utility::vector1< Size > const & r1_connids( res1.connections_to_residue( res2 ) );
		for ( Size ii = 1; ii <= r1_connids.size(); ++ii ) {
			++n_connect_;

			Size const r1_conn_id = r1_connids[ ii ];
			Size const r2_conn_id = res1.actual_residue_connection( r1_conn_id ).connid();

			Size const r1conn_atom = res1.residue_connection( r1_conn_id ).atomno();
			Size const r2conn_atom = res2.residue_connection( r2_conn_id ).atomno();

			res1_conn_point_path_dists_.push_back( & (res1.path_distance( r1conn_atom )) );
			res2_conn_point_path_dists_.push_back( & (res2.path_distance( r2conn_atom )) );
		}
	}

	if ( res1.is_pseudo_bonded( res2 ) ) {
		PseudoBondCollectionCOP pbc = res1.get_pseudobonds_to_residue( res2.seqpos() );
		for ( PseudoBondCollection::PBIter pb_iter = pbc->iter_begin(),
				pb_iter_end = pbc->iter_end(); pb_iter != pb_iter_end; ++pb_iter ) {
			++n_pconnect_;

			Size const r1_conn_id = res1.seqpos() < res2.seqpos() ? pb_iter->lr_conn_id() : pb_iter->ur_conn_id();
			Size const r2_conn_id = res1.seqpos() < res2.seqpos() ? pb_iter->ur_conn_id() : pb_iter->lr_conn_id();

			Size const r1conn_atom = res1.residue_connection( r1_conn_id ).atomno();
			Size const r2conn_atom = res2.residue_connection( r2_conn_id ).atomno();

			res1_pbconn_point_path_dists_.push_back( & (res1.path_distance( r1conn_atom )) );
			res2_pbconn_point_path_dists_.push_back( & (res2.path_distance( r2conn_atom )) );
			pb_lengths_.push_back( pb_iter->nbonds() );
		}
	}
	debug_assert( n_connect_ == res1_conn_point_path_dists_.size() );
	debug_assert( n_connect_ == res2_conn_point_path_dists_.size() );
	debug_assert( n_pconnect_ == res1_pbconn_point_path_dists_.size() );
	debug_assert( n_pconnect_ == res2_pbconn_point_path_dists_.size() );
}

/// @brief Create a count pair object that pretends there exist
/// a chemical bond between some number of atoms in residue 1 and
/// some number of atoms in residue2; the bond_pairs vector is
/// a set of ordered-pairs of atom-indices, where the first
/// is an atom from restype1 and the second is an atom of restype2.
CountPairGeneric::CountPairGeneric(
	chemical::ResidueType const & restype1,
	chemical::ResidueType const & restype2,
	utility::vector1< std::pair< Size, Size > > bond_pairs
)
:
	n_pconnect_( 0 ),
	crossover_( 3 )
{
	n_connect_ = bond_pairs.size();
	res1_conn_point_path_dists_.reserve( n_connect_ );
	res2_conn_point_path_dists_.reserve( n_connect_ );
	for ( Size ii = 1; ii <= n_connect_; ++ii ) {
		res1_conn_point_path_dists_.push_back( & (restype1.path_distance( bond_pairs[ ii ].first  )) );
		res2_conn_point_path_dists_.push_back( & (restype2.path_distance( bond_pairs[ ii ].second )) );
		//std::cout << "Bond between " << restype1.name() << " " << restype1.atom_name( bond_pairs[ ii ].first );
		//std::cout << " and " << restype2.name() << " " << restype2.atom_name( bond_pairs[ ii ].second ) << std::endl;
		//std::cout << "path distances res1: " << std::endl;
		//for ( Size jj = 1; jj <= restype1.natoms(); ++jj ) {
		// std::cout << "    " << restype1.atom_name( jj ) << " " << (*res1_conn_point_path_dists_[ ii ])[ jj ] << std::endl;
		//}
		//std::cout << "path distances res2: " << std::endl;
		//for ( Size jj = 1; jj <= restype2.natoms(); ++jj ) {
		// std::cout << "    " << restype2.atom_name( jj ) << " " << (*res2_conn_point_path_dists_[ ii ])[ jj ] << std::endl;
		//}
	}
}

CountPairGeneric::~CountPairGeneric() {}

// Generic crossover behavior must be specified
// Unneccessary use of the letter x
void
CountPairGeneric::set_crossover( Size xover )
{
	crossover_ = xover;
}

bool
CountPairGeneric::count(
	int const at1,
	int const at2,
	Real & weight,
	Size & path_dist
) const
{
	return operator()( at1, at2, weight, path_dist );
}

void
CountPairGeneric::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}

void
CountPairGeneric::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone(
		res1, res2, etable_energy, *this, emap );
}


void
CountPairGeneric::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole(
		res1, res2, etable_energy, *this, emap );
}


void
CountPairGeneric::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}


void
CountPairGeneric::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::TableLookupEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}

void
CountPairGeneric::residue_atom_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy( res1, res2, etable_energy, *this, emap );
}

void
CountPairGeneric::residue_atom_pair_energy_sidechain_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_backbone(
		res1, res2, etable_energy, *this, emap );
}


void
CountPairGeneric::residue_atom_pair_energy_sidechain_whole(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_whole(
		res1, res2, etable_energy, *this, emap );
}

void
CountPairGeneric::residue_atom_pair_energy_backbone_backbone(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_backbone_backbone( res1, res2, etable_energy, *this, emap );
}

void
CountPairGeneric::residue_atom_pair_energy_sidechain_sidechain(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	etable::AnalyticEtableEvaluator const & etable_energy,
	EnergyMap & emap
) const
{
	inline_residue_atom_pair_energy_sidechain_sidechain( res1, res2, etable_energy, *this, emap );
}


} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

