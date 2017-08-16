// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/ExtendChainMover.hh
/// @brief The ExtendChainMover Protocol
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_movers_ExtendChainMover_hh
#define INCLUDED_protocols_denovo_design_movers_ExtendChainMover_hh

// Unit headers
#include <protocols/denovo_design/movers/ExtendChainMover.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/moves/Mover.hh>

// Core headers

// C++ headers

namespace protocols {
namespace denovo_design {
namespace movers {

class ExtendChainMover : public protocols::moves::Mover {
public:
	/// @brief default constructor
	ExtendChainMover();

	/// @brief virtual constructor to allow derivation
	virtual ~ExtendChainMover();


	// mover virtuals
public:
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Parses the ExtendChainTags
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;


	protocols::moves::MoverOP
	fresh_instance() const override;

	protocols::moves::MoverOP
	clone() const override;

	// public methods
public:
	architects::DeNovoArchitect const &
	architect() const;

	void
	set_dry_run( bool const dry_run );

	void
	set_prepend( bool const prepend );

	void
	set_chain( core::Size const chain );

	void
	set_segment_names( std::string const & segment_names_str );

	void
	set_segment_names( SegmentNames const & seg_names );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	/// @brief designs the extension
	architects::DeNovoArchitectOP architect_;
	SegmentNames segment_names_;
	core::Size chain_;
	bool prepend_;
	bool dry_run_;
};

} // movers
} // denovo_design
} // protocols

#endif
