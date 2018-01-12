// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/ncbb/hbs/HbsRandomSmallMover.hh
/// @brief
/// @author Andrew Watkins, amw579@nyu.edu

#ifdef NOT_IN_SCONS_DEPRECATED

#ifndef INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsRandomSmallMover_hh
#define INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsRandomSmallMover_hh
// Unit Headers
#include <protocols/ncbb/a3b_hbs/A3BHbsRandomSmallMover.fwd.hh>
#include <protocols/ncbb/a3b_hbs/A3BHbsMover.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>
//#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {
namespace a3b_hbs {


/// @details
class A3BHbsRandomSmallMover : public protocols::moves::Mover {

public:

	/// @brief
	A3BHbsRandomSmallMover(
	);

	//HbsRandomSmallMover( core::Size hbs_position, core::Size hbs_length);
	A3BHbsRandomSmallMover( core::Size hbs_position);//, core::Size hbs_length, core::Real max_small_angle );
	A3BHbsRandomSmallMover( core::Size hbs_position, core::Real max_small_angle );

	virtual ~A3BHbsRandomSmallMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual void set_max_small_angle( core::Real angle ) { max_small_angle_ = angle; }

private:

	core::Size const hbs_seq_position_;
	core::Size const hbs_length_;
	core::Real max_small_angle_;

};//end HbsRandomSmallMover


}//namespace hbs
}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_hbs_HbsRandomSmallMover_hh

#endif
