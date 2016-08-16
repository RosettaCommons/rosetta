// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/triazolamer/TriazolamerMover.hh
/// @brief
/// @author Kevin Drew, kdrew@nyu.edu
#ifndef INCLUDED_protocols_simple_moves_triazolamer_TriazolamerMover_hh
#define INCLUDED_protocols_simple_moves_triazolamer_TriazolamerMover_hh
// Unit Headers
#include <protocols/simple_moves/triazolamer/TriazolamerMover.fwd.hh>
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
class TriazolamerMover : public protocols::moves::Mover {

public:

	/// @brief
	TriazolamerMover( core::Size triazolamer_seq_position );
	TriazolamerMover( core::Size triazolamer_seq_position, core::Real phi_angle, core::Real psi_angle );

	virtual ~TriazolamerMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual void set_phi( core::Real angle ) { phi_angle_ = angle; }
	virtual void set_psi( core::Real angle ) { psi_angle_ = angle; }

	virtual void update_hydrogens( core::pose::Pose & pose ) { update_hydrogens_( pose ); }

private:

	core::Size const triazolamer_pre_pos_;
	core::Size const triazolamer_post_pos_;
	core::Real phi_angle_;
	core::Real psi_angle_;

private:

	void update_hydrogens_( core::pose::Pose & pose );

};//end TriazolamerMover


}//namespace triazolamer
}//namespace simple_moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_triazolamer_TriazolamerMover_hh
