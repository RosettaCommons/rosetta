// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief  swaps both fragment xyz and amino acid
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_moves_ResTypeFragmentMover_HH
#define INCLUDED_protocols_simple_moves_ResTypeFragmentMover_HH

// Unit Headers
#include <protocols/simple_moves/ResTypeFragmentMover.fwd.hh>

// Package Headers
#include <protocols/simple_moves/FragmentMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class ResTypeFragmentMover : virtual public ClassicFragmentMover {
public:
	typedef ClassicFragmentMover Parent;

public:
	ResTypeFragmentMover(
		core::fragment::FragSetCOP fragset);


	ResTypeFragmentMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap);

	~ResTypeFragmentMover();

	ResTypeFragmentMoverOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<ResTypeFragmentMover>( ClassicFragmentMover::shared_from_this() ); }

protected:

	// frame and fragment of choice, returns false if no good fragment is found
	virtual bool apply_frames( core::pose::Pose &pose, core::fragment::FrameList const& frames ) const ;
	void swap_residue_types( core::pose::Pose &pose, std::string const & sequence, core::Size const startSeqPos ) const ;
};

} //simple_moves
} //protocols

#endif
