// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/GlycanInfoMover.hh
/// @brief Simple class for outputting glycan information. Currently, it simply prints the information.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_analysis_GlycanInfoMover_HH
#define INCLUDED_protocols_analysis_GlycanInfoMover_HH

// Unit headers
#include <protocols/analysis/GlycanInfoMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace analysis {

///@brief Simple class for outputting glycan information. Currently, it simply prints the information.
class GlycanInfoMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	GlycanInfoMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	GlycanInfoMover( GlycanInfoMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~GlycanInfoMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;
	
	void
	apply_const (core::pose::Pose const & pose );
	
	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	bool
	reinitialize_for_each_job() const override {
		return true;
	}
	
public:
	
	///@brief Get a string of where the attachment is occuring (1/4 etc.)
	std::string
	get_attachment_point_string( core::pose::Pose const & pose, core::Size resnum);
	
	

	//GlycanInfoMover & operator=( GlycanInfoMover const & src );

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

private: // methods

private: // data

};

std::ostream &
operator<<( std::ostream & os, GlycanInfoMover const & mover );

} //protocols
} //analysis

#endif //protocols_analysis_GlycanInfoMover_HH
