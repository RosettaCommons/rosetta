// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleLibrary.hh
/// @brief  Molecular mechanics bond angle library class
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleLibrary_hh
#define INCLUDED_core_scoring_mm_MMBondAngleLibrary_hh

// Unit headers
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.fwd.hh>


namespace core {
namespace scoring {
namespace mm {

// all ints for now
typedef utility::keys::Key3Tuple< int, int, int > mm_bondangle_atom_tri;
typedef utility::keys::Key2Tuple< double, double > mm_bondangle_param_set;
typedef std::multimap< mm_bondangle_atom_tri, mm_bondangle_param_set > mm_bondangle_library;
typedef std::multimap< mm_bondangle_atom_tri, mm_bondangle_param_set >::const_iterator mm_bondangle_library_citer;
typedef std::pair< mm_bondangle_library_citer, mm_bondangle_library_citer > mm_bondangle_library_citer_pair;

class MMBondAngleLibrary  : public utility::pointer::ReferenceCount
{

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMBondAngleLibrary();
	/// @brief ctor
	MMBondAngleLibrary( std::string filename, core::chemical::MMAtomTypeSetCAP mm_atom_set );

	/// @brief lookup by atom type ints
	mm_bondangle_library_citer_pair
	lookup( int atom1, int atom2, int atom3 ) const;

	/// @brief lookup by atom type strings
	mm_bondangle_library_citer_pair
	lookup( std::string atom1, std::string atom2,
		std::string atom3 ) const;

	/// @brief pretty print / debug
	void pretty_print() const;
	void pretty_print( int atom1, int atom2, int atom3 ) const;
	void pretty_print(  std::string atom1, std::string atom2,
		std::string atom3 ) const;

private:

	/// @brief library that contains bond angle params for sets in which all mm atom types
	/// are specified for all 3 positions
	mm_bondangle_library fully_assigned_mm_bondangle_library_;

	/// @brief library that contains bond angle params for sets in which all mm atom types
	/// are NOT specified for all 3 positions
	mm_bondangle_library wildcard_mm_bondangle_library_;

	core::chemical::MMAtomTypeSetCAP mm_atom_set_;

};


} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_bondangle_library_HH
