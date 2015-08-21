// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /devel/LoopExtend/LoopExtendMover.hh
/// @brief protocol-level mover for loop extension
/// @author Daniel J. Mandell

#ifndef INCLUDED_devel_loop_extend_LoopExtendMover_hh
#define INCLUDED_devel_loop_extend_LoopExtendMover_hh

// Unit Headers
#include <devel/loop_extend/LoopExtendMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <core/types.hh>

#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


namespace devel {
namespace loop_extend {

class LoopExtendMover : public protocols::moves::Mover {

public:
	//@brief constructor with arguments

	//virtual ~LoopExtendMover();
	//LoopExtendMover() {}

	LoopExtendMover(protocols::loops::Loop loop, core::Size extend_len) {
		loop_ = loop;
		extend_len_ = extend_len;
		init_phi_ = -150.0;
		init_psi_ = 150.0;
		init_omega_ = 180.0;
	}

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:

	protocols::loops::Loop loop_;
	core::Size extend_len_;
	core::Real init_phi_;
	core::Real init_psi_;
	core::Real init_omega_;

}; //class LoopExtendMover


}//LoopExtend
}//devel

#endif //INCLUDED_devel_LoopExtend_LoopExtendMover_FWD_HH
