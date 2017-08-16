// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/EnvMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_EnvMover_hh
#define INCLUDED_protocols_environment_EnvMover_hh

// Unit Headers
#include <protocols/environment/EnvMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/environment/Environment.fwd.hh>
#include <protocols/environment/ClientMover.hh>

// Project headers
#include <protocols/moves/MoverContainer.hh>

// C++ Headers
#include <set>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EnvMover : public moves::Mover {

public:
	EnvMover();

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const& pose ) override;

	~EnvMover() override;

	void apply( Pose& pose ) override;

	void add_apply_mover( protocols::moves::MoverOP );

	void add_registered_mover( protocols::moves::MoverOP );

	Environment& env() { return *env_; }



	moves::MoverOP clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	void parse_subtag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const& pose );

	EnvironmentOP env_;
	moves::SequenceMoverOP movers_;
	std::set< moves::MoverOP > reg_only_movers_;

}; // end EnvMover base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvMover_hh
