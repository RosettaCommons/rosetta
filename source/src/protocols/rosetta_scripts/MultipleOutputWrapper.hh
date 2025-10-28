// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/MultipleOutputWrapper.hh
/// @brief  This mover wraps another mover or a ROSETTASCRIPTS block to obtain additional,
//    related poses. A new instance of the mover or protocol is created for each iteration.
/// @author Luki Goldschmidt (lugo@uw.edu)

#ifndef INCLUDED_protocols_rosetta_scripts_MultipleOutputWrapper_hh
#define INCLUDED_protocols_rosetta_scripts_MultipleOutputWrapper_hh

// Package headers
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace rosetta_scripts {

class MultipleOutputWrapper;

class MultipleOutputWrapper : public protocols::moves::Mover {

public:
	/// @brief No-argument constructor
	MultipleOutputWrapper();

	/// @brief Virtual destructor
	~MultipleOutputWrapper() override = default;

	protocols::moves::MoverOP clone() const override {
		return utility::pointer::make_shared< MultipleOutputWrapper >(*this);
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return utility::pointer::make_shared< MultipleOutputWrapper >();
	}

	void apply(core::pose::Pose& pose) override;
	core::pose::PoseOP get_additional_output() override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string name_;
	utility::tag::TagCOP rosetta_scripts_tag_;
	core::pose::PoseOP reference_pose_;
	core::Size max_poses_;
	core::Size max_attempts_;
	core::Size n_poses_;
	bool keep_mover_state_;
	protocols::moves::MoverOP mover_;

protected:
	bool generate_pose(core::pose::Pose &);
};

} //rosetta_scripts
} //protocols

#endif
