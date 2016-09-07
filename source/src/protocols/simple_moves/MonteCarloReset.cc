// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/MonteCarloReset.cc
/// @author Sarel Fleishman (sarelf@uw.edu)


// Unit Headers
#include <protocols/simple_moves/MonteCarloReset.hh>
#include <protocols/simple_moves/MonteCarloResetCreator.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.MonteCarloReset" );

using namespace core;

namespace protocols {
namespace simple_moves {

std::string
MonteCarloResetCreator::keyname() const
{
	return MonteCarloResetCreator::mover_name();
}

protocols::moves::MoverOP
MonteCarloResetCreator::create_mover() const {
	return protocols::moves::MoverOP( new MonteCarloReset );
}

std::string
MonteCarloResetCreator::mover_name()
{
	return "MonteCarloReset";
}

std::string
MonteCarloReset::get_name() const {
	return MonteCarloResetCreator::mover_name();
}


/// @brief default constructor
MonteCarloReset::MonteCarloReset():
	Mover("MonteCarloReset"),
	MC_mover_( /* NULL */ )
{
}

/// @brief destructor
MonteCarloReset::~MonteCarloReset()= default;

/// @brief clone this object
MoverOP
MonteCarloReset::clone() const
{
	return MoverOP( new MonteCarloReset( *this ) );
}

/// @brief create this type of object
MoverOP
MonteCarloReset::fresh_instance() const
{
	return MoverOP( new MonteCarloReset() );
}

GenericMonteCarloMoverOP
MonteCarloReset::get_MC() const{
	return( MC_mover_ );
}

void
MonteCarloReset::set_MC( GenericMonteCarloMoverOP mc ){
	MC_mover_ = mc;
}

void
MonteCarloReset::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &movers, Pose const & pose ){
	std::string const mc_name( tag->getOption< std::string >( "MC_name" ) );
	auto find_mover( movers.find( mc_name ) );
	if ( find_mover == movers.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "MC mover not found by MonteCarloReset" );
	}

	set_MC( utility::pointer::dynamic_pointer_cast< protocols::simple_moves::GenericMonteCarloMover > ( find_mover->second ) );
	Pose temp_pose( pose );
	get_MC()->initialize();
	get_MC()->reset( temp_pose );
	TR<<"Setting MonteCarlo container with mover "<<mc_name<<std::endl;
}

void
MonteCarloReset::apply( core::pose::Pose & pose ){
	MC_mover_->reset( pose );
	TR<<"MC mover reset"<<std::endl;
}

} // ns simple_moves
} // ns protocols
