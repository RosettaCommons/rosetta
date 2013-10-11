// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignEndsMover.hh

#ifndef INCLUDED_devel_splice_AlignEndsMover_hh
#define INCLUDED_devel_splice_AlignEndsMover_hh

#include <devel/splice/AlignEndsMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <algorithm>
#include <core/kinematics/Jump.hh>

// C++ Headers
namespace devel {
namespace splice {

class AlignEndsMover : public protocols::moves::Mover {
public:
	AlignEndsMover();
    ~AlignEndsMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

		core::Real distance_threshold() const{ return distance_threshold_; }
		void distance_threshold( core::Real const r ){ distance_threshold_ = r;}

		core::Size neighbors() const{ return neighbors_; }
		void neighbors( core::Size n ){ neighbors_ = n; }

		core::Size N_terminal_count() const{ return N_terminal_count_; }
		void N_terminal_count( core::Size n ){ N_terminal_count_ = n; }

		bool odd() const{ return odd_; }
		void odd( bool const b ) { odd_ = b; }

		bool even() const{ return even_; }
		void even( bool const b ) { even_ = b; }

		core::pose::PoseOP template_pose() const;
		void template_pose( core::pose::PoseOP p );

		core::Size stagger() const{ return stagger_; }
		void stagger( core::Size const s ){ stagger_ = s; }
private:
  utility::vector1< core::Size > reference_positions( core::pose::Pose const & p ) const;
	core::Real distance_threshold_; // dflt 16;
	core::Size neighbors_, N_terminal_count_; //dflt 6, 3
	bool odd_, even_; // dflt true, true
	core::pose::PoseOP template_pose_; //dflt NULL
	core::Size stagger_; // dflt 0;
};


} // simple_moves
} // protocols

#endif
