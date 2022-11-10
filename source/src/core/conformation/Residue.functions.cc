// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/Residue.functions.cc
/// @brief  Implementation of non-member functions that operate on Residue objects
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/conformation/Residue.functions.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueConnection.hh>

// Project headers
#include <core/id/PartialAtomID.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

// C++ headers
#include <set>


namespace core {
namespace conformation {


void
idealize_hydrogens(
	Residue & res,
	Conformation const & conf
)
{
	//debug_assert( res.seqpos() != -1 ); // commented out: seqpos is Size, so can number be negative
	utility::vector1< bool > atoms_placed( res.natoms(), false );
	for ( Size ii = 1; ii <= res.nheavyatoms(); ++ii ) atoms_placed[ ii ] = true;
	int n_remaining = res.natoms() - res.nheavyatoms();

	/// Some hydrogen atom locations are calculated as a function of
	/// other hydrogen atom locations.  2HB on phe depends on 1HB, e.g.
	/// There must be some hydrogen whose coordinates are determined
	/// entirely by heavy atoms.  Usually, this hydrogen appears first
	/// in a list of hydrogens, but this isn't guaranteed by the ResidueType
	/// atom ordering, so this function makes several passes over the set of
	/// hydrogens (at most 3 passes assuming we never work with methane)

	while ( n_remaining > 0 ) {
		Size n_placed_this_iteration = 0;
		for ( Size ii = res.nheavyatoms() + 1; ii <= res.natoms(); ++ii ) {
			bool all_stubs_placed = true;
			for ( Size jj = 1; jj <= 3; ++jj ) {
				// assume that all atoms on neighboring residues have been placed
				// but check for all internal atoms whether they have been placed.
				if ( res.icoor( ii ).stub_atom( jj ).is_internal() &&
						! atoms_placed[ res.icoor( ii ).stub_atom( jj ).atomno() ] ) {
					all_stubs_placed = false;
					break;
				}
			}
			if ( all_stubs_placed ) {
				Vector iixyz = res.icoor( ii ).build( res, conf );
				res.set_xyz( ii, iixyz );
				++n_placed_this_iteration;
				atoms_placed[ ii ] = true;
				--n_remaining;
			}
		}
		if ( n_placed_this_iteration == 0 ) {
			std::cerr << "Error from core::conformation::Residue.functions.cc.";
			std::cerr << "Could not place the following hydrogens: ";
			for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
				if ( ! atoms_placed[ ii ] ) {
					std::cerr << res.atom_name( ii ) << " with atom stubs: ";
					std::cerr << res.atom_name( res.icoor( ii ).stub_atom(1).atomno() ) << ", ";
					std::cerr << res.atom_name( res.icoor( ii ).stub_atom(2).atomno() ) << ", & ";
					std::cerr << res.atom_name( res.icoor( ii ).stub_atom(3).atomno() ) << std::endl;
				}
			}
			utility_exit_with_message( "Failed to place ideal hydrogen positions" );
		}
	}

}

/// @details hokey "update chi from coordinates" useful for when
/// the coordinates are known for a rotamer (specifically, a residue
/// living outside of a conformation object); after the coordinates are
/// set, the chis have to be updated -- that won't happen automatically.
void set_chi_according_to_coordinates(
	conformation::Residue & rotamer
)
{
	for ( Size ii = 1; ii <= rotamer.nchi(); ++ii ) {
		chemical::AtomIndices const & ii_chi_atoms( rotamer.chi_atoms( ii ) );
		Real const ii_chi = numeric::dihedral(
			rotamer.xyz( ii_chi_atoms[ 1 ]),
			rotamer.xyz( ii_chi_atoms[ 2 ]),
			rotamer.xyz( ii_chi_atoms[ 3 ]),
			rotamer.xyz( ii_chi_atoms[ 4 ]) );
		rotamer.chi()[ ii ] = ii_chi;
	}
	if ( rotamer.is_protein() ) {}
}

/// @details Inserts the partial atom ids for atoms along the mainchain,
/// retrieving the residue connection indices for the atoms that extend
/// off into the lower- and upper residues from the ResConnIDs stored
/// in the Residue's connect_map. Note that if torsion 1 is requested,
/// then the input Residue must have a lower connection and that if a
/// torsion within 2 chemical bonds of the end of the mainchain atoms
/// for a residue is requested (e.g. psi for canonical AAs), then
/// the input Residue must have an upper connection.
///
/// @note The typical context in which this is used is when a term
/// is inserting the partial atom IDs for the mainchain torsions and
/// and these torsions are defined by overlapping sets of atoms.
/// Since the terms that are reporting these partial atom IDs should
/// probably not be reporting a single atom multiple times, inserting
/// these atoms into a std::set will produce non-redundancy.
void
insert_partial_atom_ids_for_mainchain_torsion(
	Residue const & rsd,
	Size const mainchain_torsion,
	std::set< id::PartialAtomID > & atoms
)
{
	debug_assert( mainchain_torsion > 0 );
	int mc_tor_int = static_cast<int>(mainchain_torsion);
	for ( int ii = mc_tor_int - 1; ii <= mc_tor_int + 2; ++ii ) {
		if ( ii < 1 ) {
			Size const lower_connect_ind = rsd.lower_connect().index();
			if ( rsd.connection_incomplete( lower_connect_ind ) ) continue;
			Size const conn_id_on_lower_residue = rsd.residue_connection_conn_id( lower_connect_ind );
			Size const lower_connect_rsd = rsd.residue_connection_partner( lower_connect_ind );
			atoms.insert( id::PartialAtomID(
				conn_id_on_lower_residue,
				lower_connect_rsd,
				0 ));
		} else if ( ii <= (int) rsd.mainchain_atoms().size() ) {
			atoms.insert( id::PartialAtomID( ii, rsd.seqpos() ));
		} else {
			Size const upper_connect_ind = rsd.upper_connect().index();
			if ( rsd.connection_incomplete( upper_connect_ind ) ) continue;
			Size const conn_id_on_upper_residue = rsd.residue_connection_conn_id( upper_connect_ind );
			Size const upper_connect_rsd = rsd.residue_connection_partner( upper_connect_ind );
			atoms.insert( id::PartialAtomID(
				conn_id_on_upper_residue,
				upper_connect_rsd,
				ii - 1 - (int) rsd.mainchain_atoms().size() ) );
		}
	}
}


}
}
