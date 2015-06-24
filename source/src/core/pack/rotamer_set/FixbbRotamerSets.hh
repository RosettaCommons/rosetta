// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/RotamerSet/FixbbRotamerSets.hh
/// @brief  Fixed-backbone Residue Sets interface class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_FixbbRotamerSets_hh
#define INCLUDED_core_pack_rotamer_set_FixbbRotamerSets_hh

// Unit Headers
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>

// Package Headers

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>


// Project Headers
#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <core/types.hh>

#ifdef WIN32
#include <core/pack/rotamer_set/RotamerSet.hh>
#endif

namespace core {
namespace pack {
namespace rotamer_set {


class FixbbRotamerSets : public RotamerSetsBase
{
public:
	typedef RotamerSetsBase parent;
	typedef utility::vector1< RotamerSetOP > RotamerSetVector;

public:
	FixbbRotamerSets();
	virtual ~FixbbRotamerSets();

	virtual
	RotamerSetCOP
	rotamer_set_for_residue( uint resid ) const = 0;

	virtual
	RotamerSetOP
	rotamer_set_for_residue( uint resid ) = 0;

	virtual
	RotamerSetCOP
	rotamer_set_for_moltenresidue( uint moltenresid ) const = 0;

	virtual
	RotamerSetOP
	rotamer_set_for_moltenresidue( uint moltenresid ) = 0;

	virtual
	RotamerSetVector::const_iterator begin() = 0;

	virtual
	RotamerSetVector::const_iterator end() = 0;

	virtual
	utility::vector1< uint > const &
	resid_2_moltenres_vector() const = 0;

	virtual
	utility::vector1< uint > const &
	moltenres_2_resid_vector() const = 0;

	virtual
	void
	show( std::ostream & out ) const = 0;
};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_RotamerSet_RotamerSets_HH
