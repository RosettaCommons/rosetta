// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/conformation/Residue.functions.cc
/// @brief  Implementation of non-member functions that operate on Residue objects
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/conformation/Residue.functions.hh>

// Package headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>

// Project headers
#include <core/chemical/AtomICoor.hh>
#include <core/types.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>


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
	if ( rotamer.is_protein() ) {

	}
}

}
}
