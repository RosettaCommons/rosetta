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
	~AddSegmentDataMover() override;


public:
	// mover virtual API
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//AddSegmentDataMover & operator=( AddSegmentDataMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Sets the name of the segment which will be created
	void
	set_segment_name( std::string const & name ) { segment_name_ = name; }

	/// @brief Sets the abego of the segment to be created
	void
	set_abego( std::string const & abego ) { abego_ = abego; }

	/// @brief Sets the secondary structure of the segment to be created
	void
	set_secstruct( std::string const & secstruct ) { secstruct_ = secstruct; }

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
