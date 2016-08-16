// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoopCloser.hh
///
/// @brief Mover subclass that serves only to close the specified loop
/// @author Tim Jacobs


#ifndef INCLUDED_devel_loop_creation_LoopCloser_HH
#define INCLUDED_devel_loop_creation_LoopCloser_HH

//Unit
#include <devel/loop_creation/LoopCloser.fwd.hh>

//protocols
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>

namespace devel {
namespace loop_creation {

class LoopCloser : public protocols::moves::Mover
{
public:

	/// @brief default constructor
	LoopCloser();

	/// @brief explicit constructor
	LoopCloser(bool prevent_nonloop_modifications);

	/// @brief Was the most recent loop-closure attempt successful?
	virtual bool success() const;

	void
	loop(protocols::loops::Loop loop);

	protocols::loops::Loop
	loop() const;

	void
	prevent_nonloop_modifications(bool prevent_nonloop_modifications);

	bool
	prevent_nonloop_modifications() const;

protected:

	//Was the most recent loop-closure attempt successful?
	bool success_;

private:

	//The loop to close
	protocols::loops::Loop loop_;

	//Should this mover be allowed to make modifications outside of the loop region?
	bool prevent_nonloop_modifications_;

};

} //loop creation
} //devel

#endif
