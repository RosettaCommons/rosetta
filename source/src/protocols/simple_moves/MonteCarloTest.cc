// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/MonteCarloTest.cc
/// @author Sarel Fleishman (sarelf@uw.edu)


// Unit Headers
#include <protocols/simple_moves/MonteCarloTest.hh>
#include <protocols/simple_moves/MonteCarloTestCreator.hh>
#include <protocols/simple_moves/GenericMonteCarloMover.hh>

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


// Utility headers

//// C++ headers

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MonteCarloTest" );

using namespace core;

namespace protocols {
namespace simple_moves {

std::string
MonteCarloTestCreator::keyname() const
{
	return MonteCarloTestCreator::mover_name();
}

protocols::moves::MoverOP
MonteCarloTestCreator::create_mover() const {
	return protocols::moves::MoverOP( new MonteCarloTest );
}

std::string
MonteCarloTestCreator::mover_name()
{
	return "MonteCarloTest";
}

std::string
MonteCarloTest::get_name() const {
	return MonteCarloTestCreator::mover_name();
}


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
		throw utility::excn::EXCN_RosettaScriptsOption( "MC mover not found by MonteCarloTest" );
	}

	set_MC( utility::pointer::dynamic_pointer_cast< protocols::simple_moves::GenericMonteCarloMover > ( find_mover->second ) );
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

} // ns simple_moves
} // ns protocols
