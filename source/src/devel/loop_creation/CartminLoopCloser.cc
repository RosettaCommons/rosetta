// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CartminLoopCloser.cc
///
/// @brief
/// @author Tim Jacobs

//unit
#include <devel/loop_creation/CartminLoopCloser.hh>
#include <devel/loop_creation/CartminLoopCloserCreator.hh>

//core
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//protocols
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>
#include <protocols/minimization_packing/MinMover.hh>

//basic
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

//numeric
#include <numeric/random/random.hh>

//utility
#include <utility>
#include <utility/tag/Tag.hh>

namespace devel {
namespace loop_creation {

static basic::Tracer TR( "devel.loop_creation.CartminLoopCloser" );

//****CREATOR METHODS****//
std::string
CartminLoopCloserCreator::keyname() const
{
	return CartminLoopCloserCreator::mover_name();
}

protocols::moves::MoverOP
CartminLoopCloserCreator::create_mover() const {
	return protocols::moves::MoverOP( new CartminLoopCloser );
}

std::string
CartminLoopCloserCreator::mover_name()
{
	return "CartminLoopCloser";
}
//****END CREATOR METHODS****//


/// @brief default constructor
CartminLoopCloser::CartminLoopCloser():
	scorefxn_(nullptr),
	minimization_tolerance_( 0.01 ),
	max_chainbreak_( 0.1 )
{
	init();
}

/// @brief explicit constructor
CartminLoopCloser::CartminLoopCloser(
	core::scoring::ScoreFunctionOP scorefxn,
	core::Real minimization_tolerance,
	core::Real max_chainbreak
):
	scorefxn_(std::move(scorefxn)),
	minimization_tolerance_(minimization_tolerance),
	max_chainbreak_(max_chainbreak)
{
	init();
}

protocols::moves::MoverOP
CartminLoopCloser::clone() const {
	return( protocols::moves::MoverOP( new CartminLoopCloser( *this ) ) );
}
protocols::moves::MoverOP
CartminLoopCloser::fresh_instance() const {
	return protocols::moves::MoverOP( new CartminLoopCloser );
}

std::string
CartminLoopCloser::get_name() const {
	return "CartminLoopCloser";
}

void
CartminLoopCloser::init(){
	if ( !scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}
}

void
CartminLoopCloser::apply(
	core::pose::Pose & pose
){
	core::kinematics::FoldTree saved_ft = pose.fold_tree();
	if ( prevent_nonloop_modifications() ) {
		protocols::loops::set_single_loop_fold_tree( pose, loop() );
	}
	protocols::loops::add_single_cutpoint_variant(pose, loop());

	// setup movemap
	core::kinematics::MoveMapOP mm;
	for ( Size ii=loop().start(); ii<=loop().stop(); ++ii ) {
		mm->set_bb( ii, true );

		//don't change phi for prolines
		//  if ( pose.residue(ii).aa() == chemical::aa_pro )
		//  {
		//   mm.set( id::TorsionID( id::phi_torsion, id::BB, ii ), false );
		//  }
	}
	TR.Debug << "Minimizing residues: " << loop().start() << " " << loop().stop() << std::endl;
	TR.Debug << "Cutpoint: " << loop().cut() << std::endl;

	protocols::minimization_packing::MinMoverOP min_mover(new protocols::minimization_packing::MinMover(mm, scorefxn_, "dfpmin_armijo_nonmonotone", minimization_tolerance_, false ));
	min_mover->cartesian(true);
	min_mover->apply(pose);

	//check for closure using chainbreak score
	bool closed = check_closure(pose);

	if ( closed ) {
		success_=true;
		return;
	}

	success_=false;
	core::pose::Pose saved_pose = pose;
	pose=saved_pose;
	pose.fold_tree(saved_ft);
}

bool
CartminLoopCloser::check_closure(core::pose::Pose & pose){
	core::scoring::ScoreFunctionOP chbreak_sf( new core::scoring::ScoreFunction );
	chbreak_sf->set_weight( core::scoring::chainbreak, 20.0 );

	core::Real chainbreak_score = chbreak_sf->score(pose);
	if ( chainbreak_score < max_chainbreak_ ) {
		return true;
	}
	return false;
}

/// @brief parse tag for use in RosettaScripts
void
CartminLoopCloser::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	using namespace core;

	if ( tag->hasOption("prevent_nonloop_modification") ) {
		prevent_nonloop_modifications_ = tag->getOption< bool >("prevent_nonloop_modifications");
	}

	if ( tag->hasOption("minimization_tolerance") ) {
		minimization_tolerance_ = tag->getOption< Real >("minimization_tolerance");
	}

	if ( tag->hasOption("max_chainbreak") ) {
		max_chainbreak_ = tag->getOption< Real >("max_chainbreak");
	}
}

} //loop creation
} //devel
