// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Simple app to make a movie.
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// basic headers
#include <basic/Tracer.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/exit.hh>

// core headers
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// protocols headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>

#include <protocols/moves/PyMOLMover.hh>

// devel headers
#include <devel/init.hh>

// c++ headers
#include <string>

// tracer
static basic::Tracer TR( "apps.pilot.doug" );

// simple PeptoidMoveMover class
class PeptiodMovieMover: public protocols::moves::Mover {
public:
	// default ctor
	PeptiodMovieMover(): protocols::moves::Mover("PeptiodMovieMove"){}

	// virtual dtor
	virtual ~PeptiodMovieMover(){}

	// methods
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "PeptiodMovieMover"; }
};

typedef utility::pointer::owning_ptr< PeptiodMovieMover > PeptiodMovieMoverOP;
typedef utility::pointer::owning_ptr< PeptiodMovieMover const > PeptiodMovieMoverCOP;

int
main( int argc, char* argv[] )
{

	// init command line options
	devel::init(argc, argv);

	// create an instance of PeptoidMovieMover
	PeptiodMovieMoverOP PMP( new PeptiodMovieMover() );

	// create job distributor
	protocols::jd2::JobDistributor::get_instance()->go( PMP );

	return 0;
}

void
PeptiodMovieMover::apply( core::pose::Pose & pose )
{

	using namespace core;
	using namespace protocols;

	//moves::AddPyMOLObserver( pose, true );
	moves::PyMOLMoverOP pmm ( new moves::PyMOLMover() );
	pmm->keep_history( true );

	// create score function
	TR << "Creating ScoreFunction..." << std::endl;
	core::scoring::ScoreFunctionOP score_fxn( scoring::ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS ) );

	// create a monte carlo object for everything
	TR << "Creating MC object..." << std::endl;
	moves::MonteCarloOP main_mc( new moves::MonteCarlo( pose, *score_fxn, 0.1 ) );

	TR << "Creating purturbing movers..." << std::endl;
	// get peptide start and end positions
	const core::Size pph_start( 1 );
	const core::Size pph_end( 9 );

	// create movemap for peptide
	kinematics::MoveMapOP pert_pph_mm( new kinematics::MoveMap() );
	pert_pph_mm->set_bb_true_range(pph_start, pph_end);

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pph_small( new simple_moves::SmallMover( pert_pph_mm, 0.8, 1 ) );
	pert_pph_small->set_preserve_detailed_balance( true );
	pert_pph_small->angle_max( 'H', 0.0 );
	pert_pph_small->angle_max( 'L', 0.25 );
	pert_pph_small->angle_max( 'E', 0.25 );

	simple_moves::ShearMoverOP pert_pph_shear( new simple_moves::ShearMover( pert_pph_mm, 0.8, 1 ) );
	pert_pph_shear->set_preserve_detailed_balance( true );
	pert_pph_shear->angle_max( 'H', 0.0 );
	pert_pph_shear->angle_max( 'L', 0.25 );
	pert_pph_shear->angle_max( 'E', 0.25 );

	// create random mover
	moves::RandomMoverOP pert_pph_random( new moves::RandomMover() );
	pert_pph_random->add_mover( pert_pph_small, 1 );
	pert_pph_random->add_mover( pert_pph_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pph_repeat( new moves::RepeatMover( pert_pph_random, 1000 ) );

	TR << "Setting up RT movers..." << std::endl;

	// create a task factory and tack operations
	pack::task::TaskFactoryOP pert_tf( new pack::task::TaskFactory() );

	pack::task::operation::InitializeFromCommandlineOP pert_ifcl( new pack::task::operation::InitializeFromCommandline() );
	pert_tf->push_back( pert_ifcl );

	pack::task::operation::ReadResfileOP pert_rrop( new pack::task::operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	pack::task::operation::RestrictToRepackingOP pert_rtrp( new pack::task::operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	// create a rotamer trials mover
	protocols::minimization_packing::RotamerTrialsMoverOP pert_rt(new protocols::minimization_packing::EnergyCutRotamerTrialsMover( score_fxn, pert_tf, main_mc, 0.1 /*energycut*/ ) );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_pph_repeat );
	pert_sequence->add_mover( pmm );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, main_mc ) );


	TR << "Setting up design movers..." << std::endl;

	// create a tack factory and task operations
	pack::task::TaskFactoryOP desn_tf( new pack::task::TaskFactory() );

	pack::task::operation::InitializeFromCommandlineOP desn_ifc( new pack::task::operation::InitializeFromCommandline() );
	desn_tf->push_back( desn_ifc );

	pack::task::operation::ReadResfileOP desn_rrop( new pack::task::operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	// create a pack rotamers mover
	protocols::minimization_packing::PackRotamersMoverOP desn_pr( new protocols::minimization_packing::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( score_fxn );
	desn_pr->nloop( 1 );

	TR << "Setting up minimization movers..." << std::endl;

	// create move map for minimization
	kinematics::MoveMapOP desn_mm_sc_bb( new kinematics::MoveMap() );
	kinematics::MoveMapOP desn_mm_sc( new kinematics::MoveMap() );

	// make all the residues on the peptide we are moving minimizable
	for ( Size i = 1; i <= pose.size(); ++i ) {
		desn_mm_sc_bb->set_bb( i, true );
		desn_mm_sc_bb->set_chi( i, true );
	}

	for ( Size i = 1; i <= pose.size(); ++i ) {
		desn_mm_sc->set_chi( i, true );
	}

	// create minimization mover
	minimization_packing::MinMoverOP desn_min_sc_bb( new minimization_packing::MinMover( desn_mm_sc_bb, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
	minimization_packing::MinMoverOP desn_min_sc( new minimization_packing::MinMover( desn_mm_sc, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );

	// create a sequence mover to hold pack rotamers and minimization movers
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pr );
	desn_sequence->add_mover( desn_min_sc );
	desn_sequence->add_mover( desn_min_sc_bb );

	// create a sequence mover for the main loop
	moves::SequenceMoverOP main_sequence( new moves::SequenceMover() );
	main_sequence->add_mover( pert_sequence );
	main_sequence->add_mover( pmm );
	main_sequence->add_mover( desn_sequence );
	main_sequence->add_mover( pmm );

	// create the main repeat mover
	moves::RepeatMoverOP main_repeat( new moves::RepeatMover( main_sequence, 100 ) );

	// quick min first
	//desn_pr->apply( pose );
	//desn_min_sc->apply( pose );
	//pmm->apply( pose );
	//desn_min_sc_bb->apply( pose );
	//pmm->apply( pose );

	// run the main mover
	main_repeat->apply( pose );

}
