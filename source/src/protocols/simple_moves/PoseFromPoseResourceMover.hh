// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/PoseFromPoseResourceMover.hh
/// @brief Mover that shuttles a Pose from a PoseResource into the DataMap when its parse_my_tag method is invoked
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_PoseFromPoseResourceMover_HH
#define INCLUDED_protocols_simple_moves_PoseFromPoseResourceMover_HH

// Unit headers
#include <protocols/simple_moves/PoseFromPoseResourceMover.fwd.hh>
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

///@brief Mover that shuttles a Pose from a PoseResource into the DataMap when its parse_my_tag method is invoked
class PoseFromPoseResourceMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PoseFromPoseResourceMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PoseFromPoseResourceMover( PoseFromPoseResourceMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PoseFromPoseResourceMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief No op -- this "mover" does all it's going to ever do during its parse_my_tag method.
	void
	apply( core::pose::Pose & pose ) override;

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

	//PoseFromPoseResourceMover & operator=( PoseFromPoseResourceMover const & src );

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
operator<<( std::ostream & os, PoseFromPoseResourceMover const & mover );

} //protocols
} //simple_moves

#endif //protocols_simple_moves_PoseFromPoseResourceMover_HH
