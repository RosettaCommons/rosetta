// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/PeptideStapleMover.hh
/// @brief Definition of classes for placing peptide staples in the pose.
/// @author Jacob Corn

#ifndef INCLUDED_protocols_simple_moves_PeptideStapleMover_hh
#define INCLUDED_protocols_simple_moves_PeptideStapleMover_hh

// Unit Headers
#include <protocols/simple_moves/PeptideStapleMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace simple_moves {
class PeptideStapleMover : public protocols::moves::Mover {
public:
	PeptideStapleMover( core::Size const staple_start, core::Size const staple_gap );


	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return( protocols::moves::MoverOP( new protocols::simple_moves::PeptideStapleMover( seqpos_, staple_gap_ ) ) );
	}

private:
	// private variables
	core::Size seqpos_;
	core::Size staple_gap_;

	// private fxns
	void derive_staple_constraints_( core::pose::Pose & pose );
	void minimize_( core::pose::Pose & pose );

	//void build_staple_ ( core::pose::Pose & pose, core::Size const seqpos, core::Size const staple_gap );

}; // PeptideStapleMover

} // moves
} // protocols


#endif
