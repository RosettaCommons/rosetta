// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/BumpSelector.cc
/// @brief  bump selector class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

//Unit headers
#include <core/pack/rotamer_set/BumpSelector.hh>

namespace core {
namespace pack {
namespace rotamer_set {

BumpSelector::BumpSelector() :
	starting_rot_bumpenergy_( 1e9 ),
	max_rot_bumpenergy_( 5.0 )
{}


BumpSelector::BumpSelector( Energy max_rot_energy ) :
	starting_rot_bumpenergy_( 1e9 ),
	max_rot_bumpenergy_( max_rot_energy ),
	best_rot_bumpenergy_( starting_rot_bumpenergy_ )
{}

void
BumpSelector::set_max_rot_bumpenergy( Energy setting )
{
	max_rot_bumpenergy_ = setting;
}


/// @details the bump selector is designed to include one lowest bump energy
/// rotamer for each aa-pos if (aoif) all rotamers fail the bump filter
/// this function signals the bump selector that a new seqpos+aa+aav
/// combination is being started in the rotamer creation loop, which
/// resets the selector process to the ground state
///
void
BumpSelector::reset()
{
	best_rot_bumpenergy_ = starting_rot_bumpenergy_;
}


/// @detais the bump selector is designed to include  one lowest bump energy
/// rotamer for each aa-pos if (aoif) all rotamers fail the bump filter.
/// this function decides whether the current rotamer should be deleted, or
/// if the previous best rotamer should be deleted and the current rotamer
/// should take its place
BumpSelectorDecision
BumpSelector::iterate_bump_selector(
	Energy bumpenergy
)
{

	//delete_current_rot = false;
	//delete_previous_rot = false;

	//static bool use_input_bump_cutoff( truefalseoption("pack_bump_cutoff") );
	//static float input_bump_cutoff(
	//	use_input_bump_cutoff ? realafteroption("pack_bump_cutoff") : max_rot_bumpenergy );

	//// determine if the current rot should be erased
	////
	//// reasons:
	//// 1) it's not the first rot for this aa
	////   and this rot has an bump energy higher
	////   than both the bump threshold and the
	////   that of the best rotamer so far considered
	////
	//// ! note that this introduces an unnecessary
	////   extra rot w/ include current (?still does?)
	////
	if ( best_rot_bumpenergy_ != starting_rot_bumpenergy_ &&
		( bumpenergy > max_rot_bumpenergy_ ) &&
		( bumpenergy >= best_rot_bumpenergy_) )
	{
		//delete_current_rot = true;
		return DELETE_ROTAMER;
	}

	//// determine if the previous rot should be deleted and
	//// the current rotamer inserted in its place
	////
	//// reasons:
	////  1) the current rot is not the first rot for this aa
	////     and it has a lower bump energy than the previous
	////     rot and no rotamers have a bump energy below the
	////     bump threshold
	////  2) the current rot is not the first rot for this aa
	////     and the current rot is the first rotamer with a
	////     bump energy lower than the bump threshold
	////
	if ( best_rot_bumpenergy_ != starting_rot_bumpenergy_ &&
		bumpenergy < best_rot_bumpenergy_ &&
		(best_rot_bumpenergy_ > max_rot_bumpenergy_ ))
	{
		//delete_previous_rot = true;
		best_rot_bumpenergy_ = bumpenergy;
		return DELETE_PREVIOUS_ROTAMER;
	}

	//// update best bumpenergy
	////
	if ( bumpenergy < best_rot_bumpenergy_ ) {
		best_rot_bumpenergy_ = bumpenergy;
	}
	return KEEP_ROTAMER;

}


} // rotamer_set
} // pack
} // core

