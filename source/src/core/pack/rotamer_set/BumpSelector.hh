// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/BumpSelector.cc
/// @brief  bump selector class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_rotamer_set_BumpSelector_hh
#define INCLUDED_core_pack_rotamer_set_BumpSelector_hh

#include <core/types.hh>

namespace core {
namespace pack {
namespace rotamer_set {

enum BumpSelectorDecision
{
	KEEP_ROTAMER,
	DELETE_ROTAMER,
	DELETE_PREVIOUS_ROTAMER
};

class BumpSelector
{
public:
	BumpSelector();

	/// @brief Set the max energy and "reset" the selector in a single call
	BumpSelector( Energy max_rot_energy );

	void set_max_rot_bumpenergy( Energy setting );

	/// @brief reset bump selector best energy
	void reset();

	/// @brief run bump filter for current rot
	BumpSelectorDecision
	iterate_bump_selector(
		Energy bumpenegy
	);

private:
	Energy const starting_rot_bumpenergy_;
	Energy max_rot_bumpenergy_;
	Energy best_rot_bumpenergy_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // rotamer_set
} // pack
} // core

#endif
