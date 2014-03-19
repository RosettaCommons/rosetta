// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rosetta_scripts/MultiplePoseMover.hh
/// @brief	This mover accepts multiple poses from a previous mover,
///   performs selection using a provided pose selector, 
///   applies contained ROSETTASCRIPTS protocol (ParsedProtocol),
///   and output multiple poses to the next mover of JD2.
/// @author Luki Goldschmidt (lugo@uw.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_MultiplePoseMover_hh
#define INCLUDED_protocols_rosetta_scripts_MultiplePoseMover_hh

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace rosetta_scripts {

class MultiplePoseMover;
typedef utility::pointer::owning_ptr< MultiplePoseMover > MultiplePoseMoverOP;
typedef utility::pointer::owning_ptr< MultiplePoseMover const > MultiplePoseMoverCOP;

class MultiplePoseMover : public protocols::moves::Mover {

public:
	/// @brief No-argument constructor
	MultiplePoseMover();

	/// @brief Virtual destructor
	virtual ~MultiplePoseMover() {};

	protocols::moves::MoverOP clone() const {
		return new MultiplePoseMover(*this);
	}

	protocols::moves::MoverOP fresh_instance() const {
		return new MultiplePoseMover();
	}

	void apply(core::pose::Pose& pose);
	core::pose::PoseOP get_additional_output();
	virtual std::string get_name() const;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	///@brief Used by RosettaScripts to set the previous mover to pull poses from
	void set_previous_mover( protocols::moves::MoverOP const m ) { previous_mover_ = m; }

private:
	core::Size max_poses_;
	utility::tag::TagCOP rosetta_scripts_tag_;
	utility::tag::TagCOP selector_tag_;

	protocols::moves::MoverOP previous_mover_;
	utility::vector1 < PoseSelectorOP > selectors_;

	utility::vector1 < core::pose::PoseOP > poses_;
	utility::vector1 < core::pose::PoseOP > selected_poses_;
	core::Size selected_poses_i_;

protected:
	bool process_pose( core::pose::Pose &, utility::vector1 < core::pose::PoseOP > & );
	void select_poses( utility::vector1 < core::pose::Pose > & );

};

} //rosetta_scripts
} //protocols

#endif
