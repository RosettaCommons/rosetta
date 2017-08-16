// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/BridgeChainsMover.hh
/// @brief Creates a bridge connection between two chains using remodel
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_BridgeChainsMover_hh
#define INCLUDED_protocols_denovo_design_movers_BridgeChainsMover_hh

// Unit headers
#include <protocols/denovo_design/movers/BridgeChainsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/denovo_design/connection/ConnectionArchitect.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Creates a bridge connection between two chains using remodel
class BridgeChainsMover : public protocols::moves::Mover {

public:

	BridgeChainsMover();

	/// @brief subclasses should call this constructor
	BridgeChainsMover( std::string const & class_name );

	// copy constructor (not needed unless you need deep copies)
	//BridgeChainsMover( BridgeChainsMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~BridgeChainsMover();


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
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	static void setup_attlist_for_derived_classes( utility::tag::AttributeList & attlist );

	//BridgeChainsMover & operator=( BridgeChainsMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

public:
	connection::ConnectionArchitect const &
	architect() const;

	void
	set_id( std::string const & my_id );

	void
	set_segment1_ids( std::string const & segments );

	void
	set_segment2_ids( std::string const & segments );

	void
	set_motifs( std::string const & motifs_str, std::string const & cut_resis_str );

	core::Size
	overlap() const;

	void
	set_overlap( core::Size const overlap_val );

	core::scoring::ScoreFunction const &
	scorefxn() const;

	void
	set_scorefxn( core::scoring::ScoreFunction const & sfxn );

	bool
	dry_run() const;

	void
	set_dry_run( bool const dry_run );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	connection::ConnectionArchitectOP architect_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Size overlap_;
	core::Size iterations_;
	bool dry_run_;

};

std::ostream &
operator<<( std::ostream & os, BridgeChainsMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_BridgeChainsMover_hh
