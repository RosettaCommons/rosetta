// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/MultiplePoseMover.hh
/// @brief This mover accepts multiple poses from a previous mover,
///   performs selection using a provided pose selector,
///   applies contained ROSETTASCRIPTS protocol (ParsedProtocol),
///   and output multiple poses to the next mover of JD2.
/// @author Luki Goldschmidt (lugo@uw.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_MultiplePoseMover_hh
#define INCLUDED_protocols_rosetta_scripts_MultiplePoseMover_hh

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <string>
#include <deque>

namespace protocols {
namespace rosetta_scripts {

class MultiplePoseMover;
typedef utility::pointer::shared_ptr< MultiplePoseMover > MultiplePoseMoverOP;
typedef utility::pointer::shared_ptr< MultiplePoseMover const > MultiplePoseMoverCOP;

class MultiplePoseMover : public protocols::moves::Mover {

public:
	/// @brief No-argument constructor
	MultiplePoseMover();

	/// @brief Virtual destructor
	~MultiplePoseMover() override = default;

	protocols::moves::MoverOP clone() const override {
		return utility::pointer::make_shared< MultiplePoseMover >(*this);
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return utility::pointer::make_shared< MultiplePoseMover >();
	}

	void apply(core::pose::Pose& pose) override;
	core::pose::PoseOP get_additional_output() override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	/// @brief Used by RosettaScripts to set the previous mover to pull poses from
	void set_previous_mover( protocols::moves::MoverOP const m );

	/// @brief sets rosettascripts tag
	void set_rosetta_scripts_tag( utility::tag::TagCOP tag );

	/// @brief sets the mover to use on each structure (versus the RosettaScripts script)
	void set_main_mover( protocols::moves::MoverOP const m );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	virtual bool process_pose( core::pose::Pose &, utility::vector1 < core::pose::PoseOP > & );

private:
	bool fill_input_cache();
	core::pose::PoseOP generate_pose();
	std::deque < core::pose::PoseOP > select_poses( std::deque < core::pose::PoseOP > & poses);
	std::deque < core::pose::PoseOP > process_poses( std::deque < core::pose::PoseOP > & poses);

private:
	bool cached_;
	core::Size max_input_poses_, max_output_poses_;
	utility::tag::TagCOP rosetta_scripts_tag_;
	protocols::moves::MoverOP mover_; // Use a mover instead of a RosettaScripts tag
	utility::tag::TagCOP selector_tag_;
	protocols::moves::MoverOP previous_mover_;
	utility::vector1 < PoseSelectorOP > selectors_;

	// Pose caches and counters
	std::deque < core::pose::PoseOP > pose_input_cache_, pose_output_cache_;
	core::Size poses_input_, poses_output_;

	/// @brief in its [previous] form MultiplePoseMover [was] not assignable due to presense of protocols::filters::Filters_map
	/// defined as std::map< const std::string, ...> but compiler tries to generate assigment operator anyway
	/// I don't know if that has changed now that we removed the Filters_map
	MultiplePoseMover & operator= ( const MultiplePoseMover & ) = delete;

};

} //rosetta_scripts
} //protocols

#endif
