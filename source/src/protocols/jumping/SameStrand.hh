// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief secondary structure will hold statistics about secondary structure predictions
/// sources can be from
///      - fragments
///      - psipred files ? other stuff
///
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_jumping_SameStrand_hh
#define INCLUDED_protocols_jumping_SameStrand_hh

// Unit Headers
#include <protocols/jumping/SameStrand.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
// #Include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <core/fragment/SecondaryStructure.fwd.hh>


//// C++ headers
//#include <cstdlib>
//#include <string>
//#include <vector>

namespace protocols {
namespace jumping {

/// @brief tiny helper class that knows the relative fractions of secondary structure  L,H,E
/// @detail
/// so far these fractions can be computed from a FragSet
/// other input strategies are conceivable but not implemented, yet: eg. psipred files, a bunch of poses,
class SameStrand : public utility::pointer::ReferenceCount {
public:
	/// @brief c'stor compute fractions from fragments
	SameStrand( core::fragment::SecondaryStructureOP );

	/// @brief explicit definitions of c'stor and d'stor
	virtual ~SameStrand();
	SameStrand( SameStrand const& );


	/// @brief print current choice to stream
	void show( std::ostream &os ) const;

	/// @brief new stochastic choices for strand boundaries
	void redo() const;

	/// @brief return whether residue i and j are on the same strand
	bool eval( core::Size i, core::Size j ) const;

	/// @brief ...
	Size total_residue() const {
		return total_residue_;
	}

private:

	void compute( core::fragment::SecondaryStructure const& ss ) const;
	void do_strand_sum( core::fragment::SecondaryStructure const& ss ) const;
	void do_same_strand( ) const; //uses only strand_sum_

	/// @brief store loop/strand
	mutable ObjexxFCL::FArray2D_bool same_strand_;

	/// @brief
	mutable ObjexxFCL::FArray1D_float strand_sum_;

	/// @brief length of FArrays
	core::Size total_residue_;

	/// @brief ScondaryStructure information --- needed permanently for redo() method
	core::fragment::SecondaryStructureOP secondary_structure_;

};

/// @brief output operator
inline std::ostream & operator <<(std::ostream & os, SameStrand const & t) {
	t.show( os );
	return os;
}


} //protocols
} //jumping

#endif

