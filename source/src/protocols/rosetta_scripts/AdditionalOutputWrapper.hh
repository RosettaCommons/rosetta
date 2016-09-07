// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/AdditionalOutputWrapper.cc
/// @brief This mover accepts multiple poses from a previous mover,
///   performs selection using a provided pose selector,
///   applies contained ROSETTASCRIPTS protocol (ParsedProtocol),
///   and output multiple poses to the next mover of JD2.
/// @author Luki Goldschmidt (lugo@uw.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_AdditionalOutputWrapper_hh
#define INCLUDED_protocols_rosetta_scripts_AdditionalOutputWrapper_hh

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

class AdditionalOutputWrapper;

class AdditionalOutputWrapper : public protocols::moves::Mover {

public:
	/// @brief No-argument constructor
	AdditionalOutputWrapper();

	/// @brief Virtual destructor
	~AdditionalOutputWrapper() override = default;

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new AdditionalOutputWrapper(*this) );
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return protocols::moves::MoverOP( new AdditionalOutputWrapper() );
	}

	void apply(core::pose::Pose& pose) override;
	core::pose::PoseOP get_additional_output() override;
	std::string get_name() const override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	std::string name_;
	utility::tag::TagCOP mover_tag_;
	utility::tag::TagCOP rosetta_scripts_tag_;
	core::pose::PoseOP reference_pose_;
	core::Size max_poses_;
	core::Size n_poses_;

protected:
	void generate_pose(core::pose::Pose &);
};

} //rosetta_scripts
} //protocols

#endif
