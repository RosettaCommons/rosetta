// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   PDBReloadMover.hh
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

#ifndef INCLUDED_protocols_simple_moves_PDBReloadMover_hh
#define INCLUDED_protocols_simple_moves_PDBReloadMover_hh

// Unit headers
#include <protocols/simple_moves/PDBReloadMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace simple_moves {

class PDBReloadMover : public moves::Mover {

public:
	/// @brief
	///  empty constructor
	PDBReloadMover();

	// Undefined, commenting out to fix PyRosetta build  PDBReloadMover( core::pose::Pose const & pose );

	~PDBReloadMover() override;

	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_PDBReloadMover_HH
