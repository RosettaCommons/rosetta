// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/scoring/ScoreFunctionFactory.hh> // get_score_function
#include <core/pose/Pose.fwd.hh>
#include <basic/prof.hh>


#include <core/pose/symmetry/util.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

// APL TEMP
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.symmetry.SymMinMover" );


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
	return protocols::moves::MoverOP( new SymMinMover );
}

std::string
SymMinMoverCreator::mover_name() {
	return "SymMinMover";
}

//////////////////////////
// default constructor
// proper lightweight default constructor
SymMinMover::SymMinMover()
: protocols::simple_moves::MinMover() {}

SymMinMover::SymMinMover( std::string const & name )
: protocols::simple_moves::MinMover(name) {}

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
) :
	protocols::simple_moves::MinMover(
	movemap_in, scorefxn_in, min_type_in,
	tolerance_in, use_nb_list_in,
	deriv_check_in, deriv_check_verbose_in ) {}


void
SymMinMover::apply( pose::Pose & pose )
{
	// lazy default initialization
	core::kinematics::MoveMapOP symmetric_movemap;
	if ( ! movemap() ) symmetric_movemap = core::kinematics::MoveMapOP( new MoveMap );
	else symmetric_movemap = movemap()->clone();

	apply_dof_tasks_to_movemap(pose, *symmetric_movemap);

	core::pose::symmetry::make_symmetric_movemap( pose, *symmetric_movemap ); // we do this here since this is the first time we meet the symmetric pose

	//{ // scope APL debug:
	// std::cout << "SymMinMover free DOFs" << std::endl;
	// core::conformation::symmetry::SymmetricConformation const & symm_conf (
	//  dynamic_cast< core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()) );
	// core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	// for (int jump_nbr = 1; jump_nbr <= (int)pose.num_jump(); ++jump_nbr) {
	//  if ( symmetric_movemap->get_jump( jump_nbr ) ) {
	//   std::cout << "SymMinMover jump " << jump_nbr << " is free" << std::endl;
	//  }
	// }
	//}


	if ( ! score_function() ) score_function( get_score_function() ); // get a default (INITIALIZED!) ScoreFunction

	PROF_START( basic::MINMOVER_APPLY );
	if ( !cartesian( ) ) {
		//TR << "Before minimization" << std::endl;
		//score_function()->show( TR, pose );
		//TR << std::endl;
		core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
		(*score_function())(pose);
		minimizer.run( pose, *symmetric_movemap, *score_function(), *min_options() );
		//TR << "After minimization" << std::endl;
		//score_function()->show( TR, pose );
		//TR << std::endl;
	} else {
		core::optimization::CartesianMinimizer minimizer;
		(*score_function())(pose);
		minimizer.run( pose, *symmetric_movemap, *score_function(), *min_options() );
	}
	PROF_STOP( basic::MINMOVER_APPLY );
}

std::string
SymMinMover::get_name() const {
	return SymMinMoverCreator::mover_name();
}

protocols::moves::MoverOP SymMinMover::clone() const { return protocols::moves::MoverOP( new  SymMinMover( *this ) ); }
protocols::moves::MoverOP SymMinMover::fresh_instance() const { return protocols::moves::MoverOP( new  SymMinMover ); }

void SymMinMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
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
