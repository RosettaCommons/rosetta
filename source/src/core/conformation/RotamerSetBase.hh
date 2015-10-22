// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/RotamerSetBase.hh
/// @brief  Abstract base class for a class that holds disembodied Residues (ie not contained within a Pose) as Rotamers
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_conformation_RotamerSetBase_hh
#define INCLUDED_core_conformation_RotamerSetBase_hh

// Unit Headers
#include <core/conformation/RotamerSetBase.fwd.hh>

// Package Headers
#include <core/conformation/AbstractRotamerTrie.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#ifdef WIN32
#include <core/conformation/AbstractRotamerTrie.hh> // WIN32 INCLUDE
#include <core/conformation/Residue.hh> // WIN32 INCLUDE
#endif

#include <basic/datacache/BasicDataCache.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace conformation {

class RotamerSetBase : public utility::pointer::ReferenceCount
{
private: // typedefs
	typedef utility::pointer::ReferenceCount parent;

public: // typedefs
	typedef basic::datacache::BasicDataCache   BasicDataCache;
	typedef basic::datacache::BasicDataCacheOP BasicDataCacheOP;

public:
	RotamerSetBase();
	virtual ~RotamerSetBase();

	virtual
	Size
	get_n_residue_types() const = 0;

	virtual
	Size
	get_residue_type_begin( Size which_restype ) const = 0;

	virtual
	Size
	get_n_rotamers_for_residue_type( Size which_restype ) const = 0;


	/// @brief Rotamers i to i+j of all the same residue type are grouped together.
	/// This function returns the index of the residue type in a contiguous block
	/// of rotamers.  E.g. rotamers 100 to 120 might all be lysine rotamers, and might
	/// be the 8th residue type, with the first 7 residue types spanning rotamers 1 to 99.
	/// If new lysine rotamers are appended to the end of the rotamer set, they are
	/// considered to be in a separate residue type block.  Lysine rotamers 200 to 210 might
	/// be block 15 while lysine rotamers 100 to 120 are still block 7.
	virtual
	Size
	get_residue_type_index_for_rotamer( Size which_rotamer ) const = 0;

	virtual
	Size
	num_rotamers() const = 0;

	virtual
	Size resid() const = 0;

	virtual
	conformation::ResidueCOP
	rotamer( Size rot_id ) const = 0;

	virtual
	conformation::Residue const &
	rotamer_ref( Size rot_id ) const = 0;

	/// @brief mutatable access to a single rotamer in the set.
	virtual
	conformation::ResidueOP
	nonconst_rotamer( Size rot_id ) = 0;

	virtual
	void
	store_trie( Size method_enum_id, AbstractRotamerTrieOP trie ) = 0;

	virtual
	AbstractRotamerTrieCOP
	get_trie( Size method_enum_id ) const = 0;

	/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
	BasicDataCache &
	data();

	/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
	BasicDataCache const &
	data() const;


private:
	// deny use of the copy constructor (no pass-by-value)
	RotamerSetBase( RotamerSetBase const & );

	/// @brief BasicDataCache indexed by enum in core/pack/rotamer_set/RotamerSetCacheableDataType.hh
	/// @warning DataCache must always be initialized with the number of cacheable
	///  data types -- see the last enum entry.
	BasicDataCacheOP data_cache_;

};


} // namespace conformation
} // namespace core

#endif
