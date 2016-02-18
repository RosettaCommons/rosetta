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

// Unit Headers
#include <core/scoring/trie/CPDataCorrespondence.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/PseudoBond.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace trie {

CPDataCorrespondence::CPDataCorrespondence()
:
	n_entries_( 0 ),
	max_connpoints_for_residue_( 0 ),
	has_pseudobonds_( false )
{}

void CPDataCorrespondence::n_entries( Size nentries )
{
	n_entries_ = nentries;
	entry_2_resid_.resize( n_entries_ );
	residue_connid_for_entry_connid_.resize( n_entries_ );
	std::fill( entry_2_resid_.begin(), entry_2_resid_.end(), 0 );
	nconnections_for_entry_.resize( n_entries_ );
	std::fill( nconnections_for_entry_.begin(), nconnections_for_entry_.end(), 0 );
	residue_connid_for_entry_connid_.resize( n_entries_ );
}

void CPDataCorrespondence::resid_for_entry( Size entry, Size resid )
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	entry_2_resid_[ entry ] = resid;
}

void CPDataCorrespondence::n_connpoints_for_entry( Size entry, Size nconnections )
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	debug_assert( nconnections_for_entry_[ entry ] == 0 );// call this function at most once for input value of entry
	nconnections_for_entry_[ entry ] = nconnections;
	residue_connid_for_entry_connid_[ entry ].resize( nconnections );
	std::fill(
		residue_connid_for_entry_connid_[ entry ].begin(),
		residue_connid_for_entry_connid_[ entry ].end(),
		0
	);
	if ( max_connpoints_for_residue_ < nconnections ) {
		max_connpoints_for_residue_ = nconnections;
	}
}

void CPDataCorrespondence::connid_for_entry_connpoint( Size entry, Size connpoint, Size residue_connid )
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	debug_assert( connpoint > 0 && connpoint <= nconnections_for_entry_[ entry ] );
	residue_connid_for_entry_connid_[ entry ][ connpoint ] = residue_connid;
}


Size CPDataCorrespondence::n_entries() const
{
	return n_entries_;
}

Size CPDataCorrespondence::resid_for_entry( Size entry ) const
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	debug_assert( entry_2_resid_[ entry ] != 0 );
	return entry_2_resid_[ entry ];
}

Size CPDataCorrespondence::n_connpoints_for_entry( Size entry ) const
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	debug_assert( nconnections_for_entry_[ entry ] != 0 );
	return nconnections_for_entry_[ entry ];
}

Size CPDataCorrespondence::connid_for_entry_connpoint( Size entry, Size connpoint ) const
{
	debug_assert( entry > 0 && entry <= n_entries_ );
	debug_assert( connpoint > 0 && connpoint <= nconnections_for_entry_[ entry ] );
	return residue_connid_for_entry_connid_[ entry ][ connpoint ];
}

/// ASSUMPTION: Connection topology is consistent across all rotamers.  Any bond found
/// on rotamer 1 is present on rotamer i.  Any pseudobond on rotamer 1 is found on rotamer i.
CPDataCorrespondence
create_cpdata_correspondence_for_rotamerset(
	conformation::RotamerSetBase const & rotset
)
{

	if ( rotset.num_rotamers() == 0 ) { CPDataCorrespondence cpdat; return cpdat; }
	return create_cpdata_correspondence_for_rotamer( *rotset.rotamer( 1 ) );
}

CPDataCorrespondence
create_cpdata_correspondence_for_rotamer(
	conformation::Residue const & example_rotamer
)
{
	using namespace conformation;

	Size n_connection_partners( 0 );
	std::map< Size, Size > resid_2_connection_partner_id;
	utility::vector1< Size > connection_partner_id_2_resid;
	std::map< Size, std::map< Size, Size > > bonds_to_residues;
	////  < other-resid < lr-conn-id, this-connid > >

	/// NOTE: the psuedobond collection on pseudobonded residues i and j are
	/// the same objects; its safe to take pseudobonds in the order they appear
	std::map< Size, utility::vector1< Size > > pseudobonds_to_residues;
	///   < other-resid < this-connid > >


	Size seqpos = example_rotamer.seqpos();
	for ( Size ii = 1; ii <= example_rotamer.n_possible_residue_connections(); ++ii ) {
		if ( example_rotamer.connection_incomplete( ii ) ) continue;

		Size other_resid = example_rotamer.connected_residue_at_resconn( ii );
		if ( resid_2_connection_partner_id.find( other_resid ) == resid_2_connection_partner_id.end() ) {
			++n_connection_partners;
			resid_2_connection_partner_id[ other_resid ] = n_connection_partners;
			connection_partner_id_2_resid.push_back( other_resid );
			std::map< Size, Size > ii_connmap;

			/// CONVENTION sort residue connections by their order in the lower residue.
			for ( Size jj = ii; jj <= example_rotamer.n_possible_residue_connections(); ++jj ) {
				if ( other_resid != example_rotamer.connected_residue_at_resconn( jj ) ) continue;
				Size lower_res_connid = seqpos < other_resid ? jj : example_rotamer.connect_map( jj ).connid();
				ii_connmap[ lower_res_connid ] = jj;
			}
			bonds_to_residues[ other_resid ] = ii_connmap;
		}
	}

	std::map< Size, PseudoBondCollectionCOP > const & pb_map( example_rotamer.pseudobonds());
	for ( std::map< Size, PseudoBondCollectionCOP >::const_iterator
			pbc_iter = pb_map.begin(), pbc_iter_end = pb_map.end();
			pbc_iter != pbc_iter_end; ++pbc_iter ) {

		Size const other_resid = pbc_iter->first;
		PseudoBondCollectionCOP pbs = pbc_iter->second;

		if ( resid_2_connection_partner_id.find( other_resid ) == resid_2_connection_partner_id.end() ) {
			++n_connection_partners;
			resid_2_connection_partner_id[ other_resid ] = n_connection_partners;
			connection_partner_id_2_resid.push_back( other_resid );
		}

		utility::vector1< Size > pseudobond_connection_points( pbs->size(), 0 );
		Size count = 0;
		for ( PseudoBondCollection::PBIter pbiter = pbs->iter_begin(), pbiter_end = pbs->iter_end();
				pbiter != pbiter_end; ++pbiter ) {
			++count;
			PseudoBond pb = *pbiter;
			/// leaving the comparison below as "<=" in case intra-residue pseudobonds are detected!
			/// generally inappropriate for a trie but what-the-hey?
			pseudobond_connection_points[ count ] = seqpos <= other_resid ? pb.lr_conn_id() : pb.ur_conn_id();
		}
		pseudobonds_to_residues[ other_resid ] = pseudobond_connection_points;
	}

	/*
	CPDataCorrespondence cpdata_map;
	cpdata_map.n_entries( connections_to_residues_.size() );
	for ( Size ii = 1; ii <= n_connection_partners; ++ii ) {
	cpdata_map.resid_for_entry( ii, found_connections_other_resid[ ii ] );
	cpdata_map.n_connpoints_for_entry( ii, n_connections_for_resid[ ii ] );

	utility::vector1< Size > const & ii_conns(
	rotset.rotamer( found_connections_examplerots[ ii ] )->
	connections_to_residue( found_connections_other_resid[ ii ]  ) );

	for ( Size jj = 1; jj <= n_connections_for_resid[ ii ]; ++jj ) {
	cpdata_map.connid_for_entry_connpoint( ii, jj, ii_conns[ jj ] );
	}
	}
	return cpdata_map;
	*/

	CPDataCorrespondence cpdata_map;
	cpdata_map.n_entries( n_connection_partners );
	for ( Size ii = 1; ii <= n_connection_partners; ++ii ) {
		Size const other_resid = connection_partner_id_2_resid[ ii ];
		cpdata_map.resid_for_entry( ii, other_resid );

		Size ii_n_connpoints( 0 );
		if ( bonds_to_residues.find( other_resid ) != bonds_to_residues.end() ) {
			ii_n_connpoints = bonds_to_residues[ other_resid ].size();
		}
		if ( pseudobonds_to_residues.find( other_resid ) != pseudobonds_to_residues.end() ) {
			ii_n_connpoints += pseudobonds_to_residues[ other_resid ].size();
		}

		cpdata_map.n_connpoints_for_entry( ii, ii_n_connpoints );
		Size count_connpoints( 0 );
		/// First bonds
		if ( bonds_to_residues.find( other_resid ) != bonds_to_residues.end() ) {
			std::map< Size, Size > const & ii_connmap = bonds_to_residues[ other_resid ];
			for ( std::map< Size, Size >::const_iterator
					iter = ii_connmap.begin(), iter_end = ii_connmap.end();
					iter != iter_end; ++iter ) {
				cpdata_map.connid_for_entry_connpoint( ii, ++count_connpoints, iter->second );
			}
		}
		/// Second pseudobonds.
		if ( pseudobonds_to_residues.find( other_resid ) != pseudobonds_to_residues.end() ) {
			utility::vector1< Size > const & ii_pbmap = pseudobonds_to_residues[ other_resid ];
			for ( utility::vector1< Size >::const_iterator
					iter = ii_pbmap.begin(), iter_end = ii_pbmap.end();
					iter != iter_end; ++iter ) {
				cpdata_map.connid_for_entry_connpoint( ii, ++count_connpoints, *iter );
			}
			cpdata_map.note_has_pseudobonds();
		}

	}
	return cpdata_map;
}

} // namespace trie
} // namespace scoring
} // namespace core

