// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopCloser.cc
///
/// @brief Mover subclass that serves only to close the specified loop
/// @author Tim Jacobs

//Unit
#include <devel/loop_creation/LoopCloser.hh>

namespace devel {
namespace loop_creation {

	LoopCloser::LoopCloser():
		prevent_nonloop_modifications_(true) {
	}
	
	LoopCloser::LoopCloser(
		bool prevent_nonloop_modifications
	):
		prevent_nonloop_modifications_(prevent_nonloop_modifications) {
	}
	
	
	/// @brief Was the most recent loop-closure attempt successful?
	bool
	LoopCloser::success() const{
		return success_;
	}
	
	void
	LoopCloser::loop(protocols::loops::Loop loop){
		loop_=loop;
	}
	
	protocols::loops::Loop
	LoopCloser::loop() const{
		return loop_;
	}
	
	void
	LoopCloser::prevent_nonloop_modifications(
		bool prevent_nonloop_modifications)
	{
		prevent_nonloop_modifications_=prevent_nonloop_modifications;
	}
	
	bool
	LoopCloser::prevent_nonloop_modifications() const{
		return prevent_nonloop_modifications_;
	}

} //loop creation
} //devel
