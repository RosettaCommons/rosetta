// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/automorphism.cc
///
/// @brief
/// @author Ian W. Davis


#include <core/chemical/automorphism.hh>
#include <core/chemical/Atom.hh>

#include <utility/vector1.hh>


namespace core {
namespace chemical {

/// @details First automorphism is the identity, [1 2 3 4 ... N]
/// The algorithm works its way through the list of atoms one at a time,
/// for each one trying to pair it with every other possible atom
/// (including itself) that isn't already paired with some earlier atom.
/// To be paired, the atoms need only be of the same chemical type;
/// for efficiency though, we also check that they have the same number of neighbors.
/// We also check that all the edges between the current atom and previously paired atoms
/// also exist between their paired counterparts (a condition for true automorphisms).
/// If we get to a point where we can't find a partner for some atom,
/// we know some earlier pairing is wrong, so we backtrack and try something else.
/// This algorithm would be more natural to represent recursively,
/// keeping state on the stack, but by "flattening" it like this we can use
/// the backtracking machinery to "resume" our search on the next call to next(),
/// thereby iterating over all the possible automorphisms.
/// The vector curr_[] thus takes the place of the stack, and used[] serves
/// to prevent two different atoms from being paired with the same partner simultaneously.
utility::vector1<Size>
AutomorphismIterator::next()
{
	Size i = 1; // one of the atoms
	// used[i] = has atom i already been partnered with someone?
	utility::vector1<bool> used(natoms_, false);
	while ( true ) {
		Size const j = curr_[i]; // i's partner atom (at the moment)
		// We can (and do) get j > natoms_ on subsequent calls to next();
		// this should skip down to the "backtracking" code in the else block.
		if ( j <= natoms_ && !used[j] && can_pair(i, j) && edges_match(i) ) {
			used[j] = true;
			++i; // look for a partner for the next i
			if ( i > natoms_ ) {
				// They've all been partnered successfully! Yay!
				utility::vector1<Size> retval = curr_; // make a copy
				++curr_[natoms_]; // try the next atom when we come back
				return retval;
			}
		} else {
			while ( true ) {
				++curr_[i]; // try the next atom
				if ( curr_[i] > natoms_ ) {
					// Oops, there is no next atom!
					// Reset curr_[] for this atom and subsequent ones;
					// pop up to previous i and increment it instead.
					for ( Size k = i; k <= natoms_; ++k ) curr_[k] = 1;
					--i;
					// I think we'd only get i < 0 if someone calls next() again
					// after it returns the empty list; i == 0 is the proper end condition.
					if ( i <= 0 ) {
						// We've reached the end; return an empty list as a signal
						return empty_list_;
					}
					used[ curr_[i] ] = false;
				} else {
					break; // there's a next atom to try, so we're done
				}
			}// end while loop for finding next permutation to try
		}// end if/else for trying to pair i and j
	}// end infinite while()
	return empty_list_; // to make compiler happy -- never get here
}

inline
bool
AutomorphismIterator::can_pair(Size i, Size j) {
	// RM: We special case HIS here, because historically HIS and HIS_D can be paired
	// but because the nitrogens don't have the same atom type the new two-restype
	// version doesn't work properly. So for HIS we special case Nhis and Ntrp so that
	// they can match. I'm not sure if this is really what we want, but it makes things work
	// like they did under the old single residue type version.
	if ( restype_.name3() == "HIS" && restype2_.name3() == "HIS" ) {
		if ( restype_.atom(i).element_type()->get_atomic_number() == 7 &&
				restype2_.atom(j).element_type()->get_atomic_number() == 7 ) {
			return true; // For HIS, nitrogens always match other nitrogens, regardless of type or connection.
		}
	}
	return restype_.atom(i).atom_type_index() == restype2_.atom(j).atom_type_index()
		&& restype_.nbrs(i).size() == restype2_.nbrs(j).size();
}

} // namespace chemical
} // namespace core
