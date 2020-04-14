// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/AcceptToBestMover.hh
/// @brief Mover equivalent of the accept_to_best command
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_relax_AcceptToBestMover_hh
#define INCLUDED_protocols_relax_AcceptToBestMover_hh

// Unit headers
#include <protocols/relax/AcceptToBestMover.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ headers
#include <string>

namespace protocols {
namespace relax {

class AcceptToBestMover : public moves::Mover {
public:
	typedef moves::Mover parent;

public:

	AcceptToBestMover();
	AcceptToBestMover( AcceptToBestMover const & );
	~AcceptToBestMover() override;

	void
	apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	protocols::moves::MoverOP fresh_instance() const override;
	protocols::moves::MoverOP clone() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;


	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_sfxn( core::scoring::ScoreFunctionOP const & sfxn ){
		sfxn_ = sfxn;
	}

private:
	core::scoring::ScoreFunctionOP sfxn_ = nullptr;
};

}
} // protocols

#endif
