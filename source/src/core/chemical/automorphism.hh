// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/automorphism.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_chemical_automorphism_hh
#define INCLUDED_core_chemical_automorphism_hh

#include <core/chemical/ResidueType.hh>

#include <utility/vector1.hh>


// Commented by inclean daemon #include <utility/vector1.hh>
// Commented by inclean daemon #include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace chemical {


class AutomorphismIterator; // fwd declaration
typedef utility::pointer::shared_ptr< AutomorphismIterator > AutomorphismIteratorOP;
typedef utility::pointer::shared_ptr< AutomorphismIterator const > AutomorphismIteratorCOP;

///@brief Enumerates the automorphisms of a residue, which are basically
/// chemical symmetries that affect RMSD calculations.
///
///@details Common automorphisms include flipping a phenyl ring (think Phe chi2)
/// or rotating a methyl group (or CF3, if you don't care about H).
/// However, they can be much more complicated than that,
/// and some cannot be imitated by rotation about a bond.
/// Examples include labeling a -CF3 clockwise vs. counterclockwise,
/// or swapping the -CF3 branches of -C(CF3)2R.
/// See the ligand of PDB 1PQC for a reasonably complex example.
///
/// Formally, a graph automorphism is an isomorphism of that graph with itself:
/// given a graph G(V,E) and a mapping M that relabels the vertices
/// according to M(v) -> v', then M is an automorphism iff
/// (M(u),M(v)) is an edge if and only if (u,v) is an edge.
/// If the vertices are "colored" (in our case, have atom types), it must also
/// be true that M(v) and v have the same color, for all v in V.
///
/// Thus you can re-label a phenyl ring
///
///   2  3          6  5               6  3
/// 1      4  or  1      4  but not  1      4
///   6  5          2  3               2  5
///
/// because in the last case, there are new bonds 6-3 and 2-5,
/// and missing bonds 6-5 and 2-3.
///
/// See also:  OpenEye's OEChem library and its OERMSD() function.
///
class AutomorphismIterator : public utility::pointer::ReferenceCount
{
public:

	/// @brief Including H will lead to many, many more automorphisms!
	AutomorphismIterator(ResidueType const & restype, bool includeH = false):
		restype_(restype),
		restype2_(restype),
		empty_list_()
	{
		natoms_ = (includeH ? restype_.natoms() : restype_.nheavyatoms());
		curr_.assign(natoms_, 1); // = [1, 1, 1, ..., 1]
	}

	/// @brief The mapping returned will be from restype to restype2
	/// Including H will lead to many, many more automorphisms!
	AutomorphismIterator(ResidueType const & restype, ResidueType const & restype2, bool includeH = false):
		restype_(restype),
		restype2_(restype2),
		empty_list_()
	{
		natoms_ = (includeH ? restype_.natoms() : restype_.nheavyatoms());
		core::Size natoms2 = 	(includeH ? restype2_.natoms() : restype2_.nheavyatoms());
		runtime_assert( natoms_ == natoms2 );
		curr_.assign(natoms_, 1); // = [1, 1, 1, ..., 1]
	}
	virtual ~AutomorphismIterator() {}

	/// @brief Returns the next automorphism for this residue type
	/// as a vector that maps "old" atom indices to "new" atom indices.
	/// Returns an empty vector when no more automorphisms remain.
	utility::vector1<Size>
	next();

private:
	/// @brief Are atoms i and j potentially interchangeable?
	/// @details We want this check to be fast but also to eliminate
	/// as many potential pairings as possible.
	/// We currently check (1) atom types and (2) number of neighbors.
	/// We could potentially also check atom types of neighbors,
	/// but that costs a set comparison or two array sorts, so we don't right now.
	inline
	bool
	can_pair(Size i, Size j);

	/// @brief Does the current mapping preserve all edges?
	/// @details That is, if (i,j) is an edge, is (curr_[i],curr_[j]) also an edge?
	/// Checks all edges (i,j) where j < i, for the supplied i.
	/// (Edges with j > i can't be checked becaues curr_[j] isn't valid yet.)
	inline
	bool
	edges_match(Size i) {
		AtomIndices const & nbrs = restype_.nbrs(i);
		for( Size idx = 1, end = nbrs.size(); idx <= end; ++idx ) {
			Size const j = nbrs[idx];
			if( j > i ) continue;
			if( !bonded2( curr_[i], curr_[j] ) ) return false;
		}
		return true;
	}

	/// @brief Are atoms i and j bonded to each other on Restype1?
	inline
	bool
	bonded(Size i, Size j) {
		return restype_.path_distance(i,j) == 1;
	}

	/// @brief Are atoms i and j bonded to each other on Restype2?
	inline
	bool
	bonded2(Size i, Size j) {
		return restype2_.path_distance(i,j) == 1;
	}


private:
	ResidueType const & restype_;
	ResidueType const & restype2_;
	Size natoms_;
	/// curr_[i] = current partner for atom i
	utility::vector1<Size> curr_;
	utility::vector1<Size> const empty_list_;

}; // AutomorphismIterator


} // namespace chemical
} // namespace core

#endif // INCLUDED_core_chemical_automorphism_HH
