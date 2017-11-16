// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/devel/denovo_design/ConnectJumps.cc
/// @brief The ConnectJumps
/// @details
/// @author Tom Linsky


// Unit Headers
#include <devel/denovo_design/ConnectJumps.hh>
#include <devel/denovo_design/ConnectJumpsCreator.hh>

// Basic Headers
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "devel.denovo_design.ConnectJumps" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

// XRW TEMP std::string
// XRW TEMP ConnectJumpsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ConnectJumps::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ConnectJumpsCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ConnectJumps() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ConnectJumps::mover_name()
// XRW TEMP {
// XRW TEMP  return "ConnectJumps";
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  ConnectJumps main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
ConnectJumps::ConnectJumps() :
	protocols::denovo_design::movers::BridgeChainsMover()
{
}

/// @brief destructor - this class has no dynamic allocation, so
ConnectJumps::~ConnectJumps() = default;


/// Return a copy of ourselves
protocols::moves::MoverOP
ConnectJumps::clone() const
{
	return protocols::moves::MoverOP( new ConnectJumps(*this) );
}

// XRW TEMP std::string
// XRW TEMP ConnectJumps::get_name() const
// XRW TEMP {
// XRW TEMP  return ConnectJumps::mover_name();
// XRW TEMP }

void ConnectJumps::apply( core::pose::Pose & pose )
{
	TR.Warning << "*****************************************************************************" << std::endl;
	TR.Warning << "WARNING: ConnectJumps is deprecated. Please use BridgeChains instead." << std::endl;
	TR.Warning << "*****************************************************************************" << std::endl;
	protocols::denovo_design::movers::BridgeChainsMover::apply( pose );
}

std::string ConnectJumps::get_name() const {
	return mover_name();
}

std::string ConnectJumps::mover_name() {
	return "ConnectJumps";
}

void ConnectJumps::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	BridgeChainsMover::setup_attlist_for_derived_classes( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "The ConnectJumps, as the Doxygen calls it, is deprecated. Use BridgeChains instead.", attlist );
}

std::string ConnectJumpsCreator::keyname() const {
	return ConnectJumps::mover_name();
}

protocols::moves::MoverOP
ConnectJumpsCreator::create_mover() const {
	return protocols::moves::MoverOP( new ConnectJumps );
}

void ConnectJumpsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ConnectJumps::provide_xml_schema( xsd );
}


} // namespace denovo_design
} // namespace devel
