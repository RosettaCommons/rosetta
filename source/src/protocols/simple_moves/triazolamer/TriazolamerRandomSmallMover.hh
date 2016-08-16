// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/triazolamer/TriazolamerRandomSmallMover.hh
/// @brief
/// @author
#ifndef INCLUDED_protocols_simple_moves_triazolamer_TriazolamerRandomSmallMover_hh
#define INCLUDED_protocols_simple_moves_triazolamer_TriazolamerRandomSmallMover_hh
// Unit Headers
#include <protocols/simple_moves/triazolamer/TriazolamerRandomSmallMover.fwd.hh>
#include <protocols/simple_moves/triazolamer/TriazolamerMover.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace triazolamer {


/// @details
class TriazolamerRandomSmallMover : public protocols::moves::Mover {

public:

	/// @brief
	TriazolamerRandomSmallMover(
	);

	TriazolamerRandomSmallMover( utility::vector1< core::Size > triazolamer_seq_positions );
	TriazolamerRandomSmallMover( utility::vector1< core::Size > triazolamer_seq_positions, core::Real max_small_angle );

	virtual ~TriazolamerRandomSmallMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual void set_max_small_angle( core::Real angle ) { max_small_angle_ = angle; }

private:

	utility::vector1< core::Size > const triazolamer_seq_positions_;
	core::Real max_small_angle_;

};//end TriazolamerRandomSmallMover


}//namespace triazolamer
}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_triazolamer_TriazolamerRandomSmallMover_hh
