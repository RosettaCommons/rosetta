// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/AddSegmentDataMover.hh
/// @brief Adds a segment to the structuredata
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_movers_AddSegmentDataMover_hh
#define INCLUDED_protocols_denovo_design_movers_AddSegmentDataMover_hh

// Unit headers
#include <protocols/denovo_design/movers/AddSegmentDataMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Adds a segment to the structuredata
class AddSegmentDataMover : public protocols::moves::Mover {

public:

	AddSegmentDataMover();

	// copy constructor (not needed unless you need deep copies)
	//AddSegmentDataMover( AddSegmentDataMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~AddSegmentDataMover();

	static std::string
	class_name();

public:
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//AddSegmentDataMover & operator=( AddSegmentDataMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

private:
	void
	create_segment( components::StructureData & perm ) const;

private:
	std::string segment_name_;
	std::string secstruct_;
	std::string abego_;

};

std::ostream &
operator<<( std::ostream & os, AddSegmentDataMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_AddSegmentDataMover_hh
