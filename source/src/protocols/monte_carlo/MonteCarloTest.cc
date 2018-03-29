// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/monte_carlo/MonteCarloTest.cc
/// @author Sarel Fleishman (sarelf@uw.edu)


// Unit Headers
#include <protocols/monte_carlo/MonteCarloTest.hh>
#include <protocols/monte_carlo/MonteCarloTestCreator.hh>
#include <protocols/monte_carlo/GenericMonteCarloMover.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


// Utility headers

//// C++ headers

static basic::Tracer TR( "protocols.simple_moves.MonteCarloTest" );

using namespace core;
using namespace protocols::moves;

namespace protocols {
namespace monte_carlo {

// XRW TEMP std::string
// XRW TEMP MonteCarloTestCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return MonteCarloTest::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MonteCarloTestCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MonteCarloTest );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MonteCarloTest::mover_name()
// XRW TEMP {
// XRW TEMP  return "MonteCarloTest";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MonteCarloTest::get_name() const {
// XRW TEMP  return MonteCarloTest::mover_name();
// XRW TEMP }


/// @brief default constructor
MonteCarloTest::MonteCarloTest():
	Mover("MonteCarloTest"),
	MC_mover_( /* NULL */ )
{
}

/// @brief destructor
MonteCarloTest::~MonteCarloTest()= default;

/// @brief clone this object
MoverOP
MonteCarloTest::clone() const
{
	return MoverOP( new MonteCarloTest( *this ) );
}

/// @brief create this type of object
MoverOP
MonteCarloTest::fresh_instance() const
{
	return MoverOP( new MonteCarloTest() );
}

GenericMonteCarloMoverOP
MonteCarloTest::get_MC() const{
	return( MC_mover_ );
}

void
MonteCarloTest::set_MC( GenericMonteCarloMoverOP mc ){
	MC_mover_ = mc;
}

void
MonteCarloTest::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &movers, Pose const & pose ){
	std::string const mc_name( tag->getOption< std::string >( "MC_name" ) );
	auto find_mover( movers.find( mc_name ) );
	if ( find_mover == movers.end() ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "MC mover not found by MonteCarloTest" );
	}

	set_MC( utility::pointer::dynamic_pointer_cast< protocols::monte_carlo::GenericMonteCarloMover > ( find_mover->second ) );
	Pose temp_pose( pose );
	get_MC()->initialize();
	get_MC()->reset( temp_pose );
	TR<<"Setting MonteCarlo container with mover "<<mc_name<<std::endl;
}

void
MonteCarloTest::apply( core::pose::Pose & pose ){
	bool const accept( MC_mover_->boltzmann( pose ) );
	TR<<"MC mover accept="<<accept<<std::endl;
}

std::string MonteCarloTest::get_name() const {
	return mover_name();
}

std::string MonteCarloTest::mover_name() {
	return "MonteCarloTest";
}

void MonteCarloTest::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"MC_name", xsct_rosetta_bool,
		"name of a previously defined GenericMonteCarloMover");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Associated with GenericMonteCarlo. Simply test the MC criterion of the specified "
		"GenericMonteCarloMover and save the current pose if accept. "
		"Useful in conjunction with MonteCarloRecover (below) if you're running a trajectory consisting "
		"of many different sorts of movers, and would like at each point to decide whether the pose has "
		"made an improvement.",
		attlist );
}

std::string MonteCarloTestCreator::keyname() const {
	return MonteCarloTest::mover_name();
}

protocols::moves::MoverOP
MonteCarloTestCreator::create_mover() const {
	return protocols::moves::MoverOP( new MonteCarloTest );
}

void MonteCarloTestCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MonteCarloTest::provide_xml_schema( xsd );
}


} // ns simple_moves
} // ns protocols
