// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/trie/trie.functions.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_trie_trie_functions_hh
#define INCLUDED_core_scoring_trie_trie_functions_hh

// Unit Headers
#include <core/scoring/trie/trie.functions.fwd.hh>

// Package Headers
#include <core/scoring/trie/CPDataCorrespondence.hh>
#include <core/scoring/trie/RotamerDescriptor.hh>
#include <core/scoring/trie/RotamerTrieBase.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace trie {

/// @details
/// Class AT must define a constructor that takes a residue const & and an atom index.
/// Class AT must also define an operator < and an operator ==
template < class AT, class CPDAT >
RotamerTrieBaseOP
create_trie(
	conformation::RotamerSetBase const & rotset,
	AT    const & /* dummy variable for type identification */,
	CPDAT const & /* dummy variable for type identification */,
	CPDataCorrespondence const & cpdata_map,
	Distance atomic_interaction_cutoff
)
{

	utility::vector1< RotamerDescriptor< AT, CPDAT > > rotamer_descriptors( rotset.num_rotamers() );

	for ( Size ii = 1; ii <= rotset.num_rotamers(); ++ii ) {
		conformation::ResidueCOP ii_rotamer( rotset.rotamer( ii ) );
		rotamer_descriptors[ ii ].natoms( ii_rotamer->natoms() );
		Size count_added_atoms = 0;
		for ( Size jj = 1; jj <= ii_rotamer->nheavyatoms(); ++jj ) {

			AT newatom( *ii_rotamer, jj );
			newatom.is_hydrogen( false );

			CPDAT cpdata;
			initialize_cpdata_for_atom( cpdata, jj, *ii_rotamer, cpdata_map );


			RotamerDescriptorAtom< AT, CPDAT > rdatom( newatom, cpdata );
			rotamer_descriptors[ ii ].atom( ++count_added_atoms, rdatom );

			for ( Size kk = ii_rotamer->attached_H_begin( jj ),
					kk_end = ii_rotamer->attached_H_end( jj );
					kk <= kk_end; ++kk ) {

				AT newhatom( *ii_rotamer, kk );
				newhatom.is_hydrogen( true );

				CPDAT cpdata;
				initialize_cpdata_for_atom( cpdata, kk, *ii_rotamer, cpdata_map );

				RotamerDescriptorAtom< AT, CPDAT > rdatom( newhatom, cpdata );
				rotamer_descriptors[ ii ].atom( ++count_added_atoms, rdatom );

			}
		}
		rotamer_descriptors[ ii ].rotamer_id( ii );
	}

	sort( rotamer_descriptors.begin(), rotamer_descriptors.end() );

	RotamerTrieBaseOP newtrie = RotamerTrieBaseOP( new RotamerTrie< AT, CPDAT >( rotamer_descriptors, atomic_interaction_cutoff ) );
	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		newtrie->set_resid_2_connection_entry( cpdata_map.resid_for_entry( ii ), ii );
	}
	return newtrie;
}

/// @details
/// Creates a rotamer trie given a single residue -- a one-residue trie.
/// Class AT must define a constructor that takes a residue const & and an atom index.
template < class AT, class CPDAT >
RotamerTrieBaseOP
create_trie(
	conformation::Residue const & res,
	AT    const & /* dummy variable for type identification */,
	CPDAT const & /* dummy variable for type identification */,
	CPDataCorrespondence const & cpdata_map,
	Distance atomic_interaction_cutoff
)
{

	utility::vector1< RotamerDescriptor< AT, CPDAT > > rotamer_descriptor( 1 );

	rotamer_descriptor[ 1 ].natoms( res.natoms() );
	Size count_added_atoms = 0;
	for ( Size jj = 1; jj <= res.nheavyatoms(); ++jj ) {

		AT newatom( res, jj );
		newatom.is_hydrogen( false );

		CPDAT cpdata;
		initialize_cpdata_for_atom( cpdata, jj, res, cpdata_map );


		RotamerDescriptorAtom< AT, CPDAT > rdatom( newatom, cpdata );
		rotamer_descriptor[ 1 ].atom( ++count_added_atoms, rdatom );

		for ( Size kk = res.attached_H_begin( jj ), kk_end = res.attached_H_end( jj );
				kk <= kk_end; ++kk ) {

			AT newhatom( res, kk );
			newhatom.is_hydrogen( true );

			CPDAT cpdata;
			initialize_cpdata_for_atom( cpdata, kk, res, cpdata_map );

			RotamerDescriptorAtom< AT, CPDAT > rdatom( newhatom, cpdata );
			rotamer_descriptor[ 1 ].atom( ++count_added_atoms, rdatom );
		}
	}
	rotamer_descriptor[ 1 ].rotamer_id( 1 );

	RotamerTrieBaseOP newtrie = RotamerTrieBaseOP( new RotamerTrie< AT, CPDAT >( rotamer_descriptor, atomic_interaction_cutoff ) );
	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		newtrie->set_resid_2_connection_entry( cpdata_map.resid_for_entry( ii ), ii );
	}
	return newtrie;
}


template < class CPDAT >
void
initialize_cpdata_for_atom(
	CPDAT & cp_data,
	Size atom_index,
	conformation::Residue const & res,
	CPDataCorrespondence const & cpdata_map
)
{
	for ( Size ii = 1; ii <= cpdata_map.n_entries(); ++ii ) {
		if ( res.is_bonded( cpdata_map.resid_for_entry( ii ) ) ) {
			for ( Size jj = 1; jj <= cpdata_map.n_connpoints_for_entry( ii ); ++jj ) {
				Size const connection_atom(
					res.residue_connect_atom_index( cpdata_map.connid_for_entry_connpoint( ii, jj ) ));
				cp_data.set_dist_to_connect_point(
					ii, jj, res.path_distance( atom_index, connection_atom ));
			}
		} else {
			for ( Size jj = 1; jj <= cpdata_map.n_connpoints_for_entry( ii ); ++jj ) {
				cp_data.set_dist_to_connect_point( ii, jj, etable::count_pair::INFINITE_SEPARATION );
			}
		}
	}
}

} // namespace trie
} // namespace scoring
} // namespace core

#endif
