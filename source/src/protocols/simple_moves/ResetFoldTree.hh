// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResetFoldTree.hh
/// @brief wipes out a fold tree making the first residue 0 and the last residue the length of the protein

#ifndef INCLUDED_protocols_simple_moves_ResetFoldTree_hh
#define INCLUDED_protocols_simple_moves_ResetFoldTree_hh

#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class ResetFoldTree : public moves::Mover {
public:
	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap & ,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};


} // simple_moves
} // protocols

#endif
