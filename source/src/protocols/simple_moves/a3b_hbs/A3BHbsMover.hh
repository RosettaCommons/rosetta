// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file /protocols/simple_moves/a3b_hbs/A3BHbsMover.hh
/// @brief
/// @author Andy Watkins, amw579@nyu.edu

#ifdef NOT_IN_SCONS_DEPRECATED

#ifndef INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsMover_hh
#define INCLUDED_protocols_simple_moves_a3b_hbs_A3BHbsMover_hh
// Unit Headers
#include <protocols/simple_moves/a3b_hbs/A3BHbsMover.fwd.hh>
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
class A3BHbsMover : public protocols::moves::Mover {

public:

	/// @brief
	A3BHbsMover( core::Size seq_position );
	A3BHbsMover( core::Size seq_position, core::Real phi_angle, core::Real psi_angle );

	virtual ~A3BHbsMover();
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual void set_phi( core::Real angle ) { phi_angle_ = angle; }
	virtual void set_psi( core::Real angle ) { psi_angle_ = angle; }

  core::Size seq_pos_;
private:

	core::Real phi_angle_;
	core::Real psi_angle_;

};//end A3BHbsMover


}//namespace hbs
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_hbs_HbsMover_hh

#endif

