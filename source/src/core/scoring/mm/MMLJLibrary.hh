// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJLibrary.hh
/// @brief  Molecular mechanics lj library class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMLJLibrary_hh
#define INCLUDED_core_scoring_mm_MMLJLibrary_hh

// Unit headers
#include <core/scoring/mm/MMLJLibrary.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/MMAtomTypeSet.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace mm {

// typedefs
typedef Size mm_lj_atom;
typedef utility::keys::Key2Tuple< Real, Real > mm_lj_param_set;

/// @brief A class to maintain a set of MM LJ paramaters.
///
/// @details blah
///
class MMLJLibrary  : public utility::pointer::ReferenceCount
{

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMLJLibrary();
  /// @brief ctor
  MMLJLibrary( core::chemical::MMAtomTypeSetCOP mm_atom_set );

  /// @brief blah
	inline
  mm_lj_param_set
  lookup_three_bond( Size atom ) const
	{ return mm_lj_three_bond_library_[ atom ]; }

  /// @brief blah
	inline
  mm_lj_param_set
  lookup_three_bond( std::string atom ) const
	{ return mm_lj_three_bond_library_[ mm_atom_set_.lock()->atom_type_index( atom ) ]; }

  /// @brief blah
	inline
	mm_lj_param_set
  lookup( Size atom ) const
	{ return mm_lj_library_[ atom ]; }

  /// @brief blah
	inline
	mm_lj_param_set
  lookup( std::string atom ) const
	{ return mm_lj_library_[ mm_atom_set_.lock()->atom_type_index( atom ) ]; }

	/// @brief blah
	core::chemical::MMAtomTypeSetCAP
	mm_atom_set() const
	{ return mm_atom_set_; }

  /// @brief blah
	Real
	nblist_dis2_cutoff_XX() const
	{ return nblist_dis2_cutoff_XX_; }

  /// @brief blah
	Real
	nblist_dis2_cutoff_XH() const
	{ return nblist_dis2_cutoff_XH_; }

  /// @brief blah
	Real
	nblist_dis2_cutoff_HH() const
	{ return nblist_dis2_cutoff_HH_; }

private:
  /// @brief library that contains lj params for sets in which atoms are seperated by 3 or bonds
  utility::vector1<mm::mm_lj_param_set> mm_lj_three_bond_library_;

  /// @brief library that contains lj params for sets in which atoms are seperated by 4 or more bonds
  utility::vector1<mm::mm_lj_param_set> mm_lj_library_;

  /// @brief the MMAtomTypeSet associated with the library
  core::chemical::MMAtomTypeSetCAP mm_atom_set_;

	/// @brief the cutoff distance at which the neighbor list will count two atoms as being neighbors
	/// given that those atom are both heavy, a heavy and a hydrogen, or a hydrogen and  hydrogen
	/// respectivly
	Real nblist_dis2_cutoff_XX_;
	Real nblist_dis2_cutoff_XH_;
	Real nblist_dis2_cutoff_HH_;
};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_lj_library_HH
