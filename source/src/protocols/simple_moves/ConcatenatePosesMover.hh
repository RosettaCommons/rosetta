// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ConcatenatePosesMover.hh
/// @brief links supplied Poses by their termini
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_simple_moves_ConcatenatePosesMover_HH
#define INCLUDED_protocols_simple_moves_ConcatenatePosesMover_HH

// Unit headers
#include <protocols/simple_moves/ConcatenatePosesMover.fwd.hh>
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

///@brief links supplied Poses by their termini
class ConcatenatePosesMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ConcatenatePosesMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ConcatenatePosesMover( ConcatenatePosesMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ConcatenatePosesMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	concatenate_poses(core::pose::Pose& pose);

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

	//ConcatenatePosesMover & operator=( ConcatenatePosesMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	std::string
	get_component_file() const;

	void
	set_component_file( std::string );

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private: // methods

private: // data

	std::string component_file_;

};

std::ostream &
operator<<( std::ostream & os, ConcatenatePosesMover const & mover );

} //simple_moves
} //protocols

#endif //protocols_simple_moves_ConcatenatePosesMover_HH
