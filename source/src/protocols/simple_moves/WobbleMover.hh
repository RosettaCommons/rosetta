// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


#ifndef INCLUDED_protocols_simple_moves_WobbleMover_HH
#define INCLUDED_protocols_simple_moves_WobbleMover_HH

// Unit Headers There needs to be a forward header for this class.
/// #include <protocols/abinitio/WobbleMover.fwd.hh>

// Package Headers
#include <protocols/simple_moves/SmoothFragmentMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <core/fragment/Frame.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class WobbleMover;
typedef utility::pointer::shared_ptr< WobbleMover > WobbleMoverOP;
typedef utility::pointer::shared_ptr< const WobbleMover > WobbleMoverCOP;


/// @brief A protocols::moves::Mover class for a classic-wobble analog: a smooth move followed by ccd closure
/// @detail a smooth fragment is chosen according to the FragmentCost Functor;
///   a cutpoint is inserted just in front of or just after the fragment
///   a loop is defined around the fragment and cutpoint to be closed with ccd:
///   a cut_Cterm insertion: ----lfff bbb----   f: fragment_res b: buffer_res -: immovable residues
///   a !cut_Cterm insertion: ---bbb fffl---
/// the number of b resiudes is controlled by buffer_length_ (default 3);
/// the move is used by apply() (inherited from FragmentMover).
/// the insertion and loop closure is implemented in the virtual method apply_fragment().
class WobbleMover : public protocols::simple_moves::SmoothFragmentMover {
public:
	typedef protocols::simple_moves::SmoothFragmentMover Parent;

public:
	WobbleMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap,
		protocols::simple_moves::FragmentCostOP cost
	);

	WobbleMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap );

	~WobbleMover() override;
	std::string get_name() const override;

	void set_defaults() override {
		buffer_length_=3;
		forward_threshold_ = 0.3;
		backward_threshold_ = 0.3;
	};

	void set_buffer_length( core::Size setting ) {
		buffer_length_ = setting;
	};

protected:
	bool
	apply_fragment(
		core::fragment::Frame const& frame,
		core::Size frag_num,
		core::kinematics::MoveMap const& movemap,
		core::pose::Pose &pose
	) const override;

	/// @brief close loop and return if successful ( deviations smaller than thresholds )
	bool ccd_closure(
		core::pose::Pose & pose,
		protocols::loops::Loops const & loops,
		core::kinematics::MoveMap const& mm
	) const;


private:
	/// @brief how many variable residues on opposite side of fragment
	core::Size buffer_length_;

	/// @brief cutoffs that ccd has to undercut to be accepted
	core::Real forward_threshold_;
	core::Real backward_threshold_;
};

} //abinitio
} //protocols

#endif
