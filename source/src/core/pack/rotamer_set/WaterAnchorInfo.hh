// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_core_pack_rotamer_set_WaterAnchorInfo_hh
#define INCLUDED_core_pack_rotamer_set_WaterAnchorInfo_hh

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++

namespace core {
namespace pack {
namespace rotamer_set {

class WaterAnchorInfo : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~WaterAnchorInfo();
	typedef chemical::AA AA;
	typedef chemical::ResidueType ResidueType;

public:

	Size
	anchor_residue() const;

	void
	anchor_residue( Size const rsd );

	bool
	attaches_to_residue_type( ResidueType const & type ) const;

	Size
	anchor_atom( ResidueType const & type ) const;

	void
	anchor_atom( std::string const & name );

	void
	aa( AA const & aa_in );

	Size
	nstep() const;

	void
	nstep( Size const nstep_in );

private:
	Size anchor_residue_;
	std::string anchor_atom_name_;
	AA aa_;
	Size nstep_;
};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif
