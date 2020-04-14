// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ResetFullModelInfoMover.cc
/// @brief Ensure synchronized full model info
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/simple_moves/ResetFullModelInfoMover.hh>
#include <protocols/simple_moves/ResetFullModelInfoMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.ResetFullModelInfoMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ResetFullModelInfoMover::ResetFullModelInfoMover():
	protocols::moves::Mover( ResetFullModelInfoMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ResetFullModelInfoMover::ResetFullModelInfoMover( ResetFullModelInfoMover const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ResetFullModelInfoMover::~ResetFullModelInfoMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
ResetFullModelInfoMover::apply( core::pose::Pose & pose ) {
	// If the pose has a FullModelInfo object this operation has obsoleted it.
	// Fix that, rather bluntly.
	core::pose::full_model_info::FullModelInfoOP full_model_info( new core::pose::full_model_info::FullModelInfo( pose ) );
	core::pose::full_model_info::set_full_model_info( pose, full_model_info );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ResetFullModelInfoMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ResetFullModelInfoMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void ResetFullModelInfoMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Resets a Pose's FullModelInfo", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ResetFullModelInfoMover::fresh_instance() const
{
	return utility::pointer::make_shared< ResetFullModelInfoMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ResetFullModelInfoMover::clone() const
{
	return utility::pointer::make_shared< ResetFullModelInfoMover >( *this );
}

std::string ResetFullModelInfoMover::get_name() const {
	return mover_name();
}

std::string ResetFullModelInfoMover::mover_name() {
	return "ResetFullModelInfoMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
ResetFullModelInfoMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< ResetFullModelInfoMover >();
}

std::string
ResetFullModelInfoMoverCreator::keyname() const
{
	return ResetFullModelInfoMover::mover_name();
}

void ResetFullModelInfoMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResetFullModelInfoMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ResetFullModelInfoMover const & mover )
{
	mover.show(os);
	return os;
}

} //simple_moves
} //protocols
