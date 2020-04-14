// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ResetFullModelInfoMover.hh
/// @brief Ensure synchronized full model info
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_ResetFullModelInfoMover_HH
#define INCLUDED_protocols_simple_moves_ResetFullModelInfoMover_HH

// Unit headers
#include <protocols/simple_moves/ResetFullModelInfoMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace simple_moves {

///@brief Ensure synchronized full model info
class ResetFullModelInfoMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ResetFullModelInfoMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ResetFullModelInfoMover( ResetFullModelInfoMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ResetFullModelInfoMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//ResetFullModelInfoMover & operator=( ResetFullModelInfoMover const & src );

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
operator<<( std::ostream & os, ResetFullModelInfoMover const & mover );

} //simple_moves
} //protocols

#endif //protocols_simple_moves_ResetFullModelInfoMover_HH
