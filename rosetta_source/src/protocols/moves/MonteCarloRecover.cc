// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/MonteCarloRecover.cc
/// @author Sarel Fleishman (sarelf@uw.edu)


// Unit Headers
#include <protocols/moves/MonteCarloRecover.hh>
#include <protocols/moves/MonteCarloRecoverCreator.hh>
#include <protocols/moves/GenericMonteCarloMover.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
// Parser headers
#include <protocols/moves/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>

// Utility headers

//// C++ headers

static basic::Tracer TR("protocols.moves.MonteCarloRecover");

using namespace core;

namespace protocols {
namespace moves {

std::string
MonteCarloRecoverCreator::keyname() const
{
	return MonteCarloRecoverCreator::mover_name();
}

protocols::moves::MoverOP
MonteCarloRecoverCreator::create_mover() const {
	return new MonteCarloRecover;
}

std::string
MonteCarloRecoverCreator::mover_name()
{
	return "MonteCarloRecover";
}

std::string
MonteCarloRecover::get_name() const {
	  return MonteCarloRecoverCreator::mover_name();
}


/// @brief default constructor
MonteCarloRecover::MonteCarloRecover():
	Mover("MonteCarloRecover"),
	recover_low_( true ),
	MC_mover_( NULL )
{
}

/// @brief destructor
MonteCarloRecover::~MonteCarloRecover(){}

/// @brief clone this object
MoverOP
MonteCarloRecover::clone() const
{
	return new MonteCarloRecover( *this );
}

/// @brief create this type of object
MoverOP
MonteCarloRecover::fresh_instance() const
{
	return new MonteCarloRecover();
}

GenericMonteCarloMoverOP
MonteCarloRecover::get_MC() const{
	return( MC_mover_ );
}

void
MonteCarloRecover::set_MC( GenericMonteCarloMoverOP mc ){
	MC_mover_ = mc;
}

void
MonteCarloRecover::parse_my_tag( TagPtr const tag, DataMap &, Filters_map const &, Movers_map const &movers, Pose const & pose ){
	std::string const mc_name( tag->getOption< std::string >( "MC_name" ) );
	Movers_map::const_iterator find_mover( movers.find( mc_name ) );
	if( find_mover == movers.end() )
		utility_exit_with_message( "MC mover not found by MonteCarloRecover" );

	set_MC( dynamic_cast< GenericMonteCarloMover * >( find_mover->second() ) );
	recover_low( tag->getOption< bool >( "recover_low", true ) );
	Pose temp_pose( pose );
	get_MC()->initialize();
	get_MC()->reset( temp_pose );
	TR<<"Setting MonteCarloRecover with mover "<<mc_name<<" and recover_low set to "<<recover_low()<<std::endl;
}

void
MonteCarloRecover::apply( core::pose::Pose & pose ){
	if( recover_low() )
		MC_mover_->recover_low( pose );
	else
		pose = *(MC_mover_->last_accepted_pose());
}

bool
MonteCarloRecover::recover_low() const{
	return( recover_low_ );
}

void
MonteCarloRecover::recover_low( bool const recover ){
	recover_low_ = recover;
}

} // ns moves
} // ns protocols
