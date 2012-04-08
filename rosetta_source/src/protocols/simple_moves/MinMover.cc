// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/MinMoverCreator.hh>

#include <protocols/moves/DataMap.hh>

// Package headers

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh> // getScoreFunction
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <basic/prof.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <iostream>
#include <protocols/rosetta_scripts/util.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;

static basic::Tracer TR("protocols.simple_moves.MinMover");

std::string
MinMoverCreator::keyname() const
{
	return MinMoverCreator::mover_name();
}

protocols::moves::MoverOP
MinMoverCreator::create_mover() const {
	return new MinMover;
}

std::string
MinMoverCreator::mover_name()
{
	return "MinMover";
}

// default constructor
// proper lightweight default constructor
MinMover::MinMover()
	: protocols::moves::Mover("MinMover"),
		movemap_(0),
		scorefxn_(0),
		min_options_(0),
		cartesian_(false)
		//		threshold_(1000000.0) // TODO: line can be deleted?
{
	min_options_ = new MinimizerOptions( "linmin", 0.01, true, false, false );
}

MinMover::MinMover( std::string const & name )
	: protocols::moves::Mover(name),
		movemap_(0),
		scorefxn_(0),
		min_options_(0),
		cartesian_(false)
		//		threshold_(1000000.0) // TODO: line can be deleted?
{
	min_options_ = new MinimizerOptions( "linmin", 0.01, true, false, false );
}

MinMover::~MinMover(){}

// constructor with arguments
MinMover::MinMover(
	MoveMapOP movemap_in,
	ScoreFunctionCOP scorefxn_in,
	std::string const & min_type_in,
	Real tolerance_in,
	bool use_nb_list_in,
	bool deriv_check_in /* = false */,
	bool deriv_check_verbose_in /* = false */
) : protocols::moves::Mover("MinMover"),
		movemap_( movemap_in ),
		scorefxn_( scorefxn_in ),
		min_options_(0),
		cartesian_(false),
		threshold_(1000000.0) // TODO: line can be deleted?
{
	min_options_ = new MinimizerOptions(
		min_type_in, tolerance_in, use_nb_list_in, deriv_check_in, deriv_check_verbose_in );
}

/// @brief allow non-const access to the internal minimizer options object
MinimizerOptionsOP
MinMover::min_options() {
	return min_options_;
}

/// @brief allow const access to the internal minimizer options object
MinimizerOptionsCOP
MinMover::min_options() const {
	return min_options_;
}

void
MinMover::movemap( MoveMapCOP movemap_in )
{
	runtime_assert( movemap_in );
	movemap_ = new MoveMap( *movemap_in );
}

MoveMapCOP
MinMover::movemap()
{
	return movemap_;
}

void
MinMover::score_function( ScoreFunctionCOP scorefxn_in )
{
	runtime_assert( scorefxn_in );
	scorefxn_ = scorefxn_in;
}

void
MinMover::score_function( ScoreFunction const & scorefxn_in )
{
	scorefxn_ = scorefxn_in.clone();
}

ScoreFunctionCOP
MinMover::score_function() const
{
	return scorefxn_;
}

void MinMover::min_type( std::string min_type_in ) { min_options_->min_type( min_type_in ); }
std::string MinMover::min_type() { return min_options_->min_type(); }


void MinMover::tolerance( Real tolerance_in ) { min_options_->minimize_tolerance( tolerance_in ); }
Real MinMover::tolerance() { return min_options_->minimize_tolerance(); }


void MinMover::nb_list( bool nb_list_in ) { min_options_->use_nblist( nb_list_in ); }
bool MinMover::nb_list() { return min_options_->use_nblist(); }


void MinMover::deriv_check( bool deriv_check_in ) { min_options_->deriv_check( deriv_check_in ); }
bool MinMover::deriv_check() { return min_options_->deriv_check(); }

void
MinMover::apply( pose::Pose & pose_ )
{
	// lazy default initialization
	if ( ! movemap_ ) movemap_ = new MoveMap;
	if ( ! scorefxn_ ) scorefxn_ = getScoreFunction(); // get a default (INITIALIZED!) ScoreFunction

	PROF_START( basic::MINMOVER_APPLY );
	if (!cartesian( )) {
		AtomTreeMinimizer minimizer;
		(*scorefxn_)(pose_);
		minimizer.run( pose_, *movemap_, *scorefxn_, *min_options_ );
	} else {
		CartesianMinimizer minimizer;
		(*scorefxn_)(pose_);
		minimizer.run( pose_, *movemap_, *scorefxn_, *min_options_ );
	}
	PROF_STOP( basic::MINMOVER_APPLY );

  // emit statistics
  scorefxn_->show(TR.Debug, pose_);
  TR.Debug << std::endl;
}

std::string
MinMover::get_name() const {
	return MinMoverCreator::mover_name();
}

protocols::moves::MoverOP MinMover::clone() const { return new protocols::simple_moves::MinMover( *this ); }
protocols::moves::MoverOP MinMover::fresh_instance() const { return new MinMover; }

void MinMover::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap & data,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose )
{
	if ( ! movemap_ ) movemap_ = new MoveMap;
	parse_opts( tag, data, filters, movers, pose );
	parse_chi_and_bb( tag );

/// parse_movemap will reset the movemap to minimize all if nothing is stated on it. The following section ensures that the protocol specifies a MoveMap before going that way
//	utility::vector1< TagPtr > const branch_tags( tag->getTags() );
//	utility::vector1< TagPtr >::const_iterator tag_it;
	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false );
}

void MinMover::parse_opts(
	TagPtr const tag,
	protocols::moves::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn", "score12" ) );
	score_function( data.get< ScoreFunction * >( "scorefxns", scorefxn_name ) );
	if ( tag->hasOption("jump") ) {
		if ( ! movemap_ ) movemap_ = new MoveMap;
		if( tag->getOption< core::Size > ( "jump" ) == 0 ) {
			movemap_->set_jump( false );
		} else {
			utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
			foreach(std::string jump, jumps){
				Size const value = std::atoi( jump.c_str() ); // convert to C string, then convert to integer, then set a Size (phew!)
				TR << "Setting min on jump " << value << std::endl;
				movemap_->set_jump( value, true );
			}
		}
	}
	min_type( tag->getOption< std::string >( "type", "dfpmin_armijo_nonmonotone" ) );
	tolerance( tag->getOption< core::Real >( "tolerance", 0.01 ) );

	cartesian( tag->getOption< core::Real >( "cartesian", false ) );

	//fpd if cartesian default to lbfgs minimization
	if ( cartesian() && !tag->hasOption("type") ) {
		min_type( "lbfgs_armijo_nonmonotone" );
	}
}

void MinMover::parse_chi_and_bb( TagPtr const tag )
{
	if ( ! movemap_ ) movemap_ = new MoveMap;
	bool const chi( tag->getOption< bool >( "chi" ) ), bb( tag->getOption< bool >( "bb" ) );
	movemap_->set_chi( chi );
	movemap_->set_bb( bb );
	TR<<"Options chi, bb: "<<chi<<", "<<bb<<std::endl;
}


} // moves
} // protocols
