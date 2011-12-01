// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/simple_moves/symmetry/SymMinMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

// Package headers

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh> // getScoreFunction
#include <core/pose/Pose.fwd.hh>
#include <basic/prof.hh>

// AUTO-REMOVED #include <utility/string_util.hh>

#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
// AUTO-REMOVED #include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.symmetry.SymMinMover");


namespace protocols {
namespace simple_moves {
namespace symmetry {

using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;

// creator
std::string
SymMinMoverCreator::keyname() const {
	return SymMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
SymMinMoverCreator::create_mover() const {
	return new SymMinMover;
}

std::string
SymMinMoverCreator::mover_name() {
	return "SymMinMover";
}

//////////////////////////
// default constructor
// proper lightweight default constructor
SymMinMover::SymMinMover()
	: MinMover() {}

SymMinMover::SymMinMover( std::string const & name )
	: MinMover(name) {}

SymMinMover::~SymMinMover(){}

// constructor with arguments
SymMinMover::SymMinMover(
	MoveMapOP movemap_in,
	ScoreFunctionCOP scorefxn_in,
	std::string const & min_type_in,
	Real tolerance_in,
	bool use_nb_list_in,
	bool deriv_check_in /* = false */,
	bool deriv_check_verbose_in /* = false */
) : MinMover(movemap_in, scorefxn_in, min_type_in,
					tolerance_in, use_nb_list_in,
					deriv_check_in, deriv_check_verbose_in ) {}



void
SymMinMover::apply( pose::Pose & pose )
{
	// lazy default initialization
	if ( ! movemap() ) symmetric_movemap_ =  new MoveMap;
	else symmetric_movemap_ = movemap()->clone();

	core::pose::symmetry::make_symmetric_movemap( pose, *symmetric_movemap_ ); // we do this here since this is the first time we meet the symmetric pose

	if ( ! score_function() ) score_function( getScoreFunction() ); // get a default (INITIALIZED!) ScoreFunction

	PROF_START( basic::MINMOVER_APPLY );
	if (!cartesian( )) {
		core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
		(*score_function())(pose);
		minimizer.run( pose, *symmetric_movemap_, *score_function(), *min_options() );
	} else {
		core::optimization::CartesianMinimizer minimizer;
		(*score_function())(pose);
		minimizer.run( pose, *symmetric_movemap_, *score_function(), *min_options() );
	}
	PROF_STOP( basic::MINMOVER_APPLY );
}

std::string
SymMinMover::get_name() const {
	return SymMinMoverCreator::mover_name();
}

protocols::moves::MoverOP SymMinMover::clone() const { return new  SymMinMover( *this ); }
protocols::moves::MoverOP SymMinMover::fresh_instance() const { return new  SymMinMover; }

void SymMinMover::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap & data,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose )
{
	MinMover::parse_my_tag( tag, data, filters, movers, pose );

	// symm-specific options 
}


} // symmetry
} // moves
} // protocols
