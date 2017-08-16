// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_build_LoopmodelWrapper_HH
#define INCLUDED_protocols_loop_build_LoopmodelWrapper_HH

// Unit headers
#include <protocols/loop_build/LoopmodelWrapper.fwd.hh>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loop_build {

class LoopmodelWrapper : public protocols::moves::Mover {

public:

	/// @copydoc protocols::moves::Mover::get_name

	/// @copydoc {parent}::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose) override;

	protocols::moves::MoverOP clone() const override;

	/// @brief Invoke the loopmodel app.
	/// @details This mover is a very thin wrapper around the command-line
	/// loopmodel app.  The purpose of this mover is exclusively to add support
	/// for the loopmodel application in rosetta scripts.  The behavior of loop
	/// modeling cannot be controlled at all from this mover and must instead be
	/// controlled using the usual command-line flags.
	void apply(core::pose::Pose & pose) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};

}
}

#endif
