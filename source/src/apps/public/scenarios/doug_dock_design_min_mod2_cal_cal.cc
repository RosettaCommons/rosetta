// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file apps/pilot/doug/doug_dock_design_min.cc
/// @brief
/// @author Doug Renfrew

// Here is a sketch of the basic flow of the program...
//
// Pertubation Phase
//   +-Monte Carlo Mover---------------------------------------+
//   | +-Random Mover--------------------------------------+ | |
//   | | +-Docking Mover-----------------------------------+ | |
//   | | | small rigid body movements between the peptide	 | | |
//   | | | and protein for conformational diversity        | | |
//   | | +-------------------------------------------------+ | |
//   | | +-Termini Modeling--------------------------------+ | |
//   | | | move termini with small/shear moves to generate | | |
//   | | | conformational diversity                        | | |
//   | | +-------------------------------------------------+ | |
//   | | +-Peptide Modeling--------------------------------+ | |
//   | | | move peptide with small/shear moves to generate | | |
//   | | | conformational diversity                        | | |
//   | | +-------------------------------------------------+ | |
//   | +-----------------------------------------------------+ |
//   | +-Rotamer Trials Mover--------------------------------+ |
//   | | quick sidechain packing to find optimal rotamers    | |
//   | | before the next cycle                               | |
//   | +-----------------------------------------------------+ |
//   +---------------------------------------------------------+
//
// Design Minimization Phase
//   +-Monte Carlo Mover---------------------------------------+
//   | +-Pack Rotamers Mover---------------------------------+ |
//   | | repack and design rotamers to explore sequence 	   | |
//   | | space  	                                           | |
//   | +-----------------------------------------------------+ |
//   | +-Minimization Mover----------------------------------+ |
//   | | energy minimize the current conformation before     | |
//   | | deciding to accept or reject the purtubations       | |
//   | +-----------------------------------------------------+ |
//   +---------------------------------------------------------+
//

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
// AUTO-REMOVED #include <protocols/loops/kinematic_closure/KinematicMover.hh>

// AUTO-REMOVED #include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/kinematic_closure/KinematicWrapper.hh>
// AUTO-REMOVED #include <protocols/loops/CCDLoopClosureMover.hh>

// AUTO-REMOVED #include <protocols/loops/loops_main.hh>

#include <protocols/rigid/RB_geometry.hh>
// AUTO-REMOVED #include <protocols/docking/DockingProtocol.hh>

// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
// AUTO-REMOVED #include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
// AUTO-REMOVED #include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
// AUTO-REMOVED #include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility headers
#include <devel/init.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <iostream>
// AUTO-REMOVED #include <fstream>
#include <string>
#include <sstream>

#include <utility/vector0.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <protocols/loops/Loops.fwd.hh>
#include <utility/excn/Exceptions.hh>

// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::toolbox;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer
static thread_local basic::Tracer TR( "apps.pilot.doug" );

// application specific options
namespace dddm {
// pert options
RealOptionKey const pert_mc_temp( "dddm::pert_mc_temp" );
RealOptionKey const pert_dock_rot_mag( "dddm::pert_dock_rot_mag" );
RealOptionKey const pert_dock_trans_mag( "dddm::pert_dock_trans_mag" );
RealOptionKey const pert_pep_small_temp( "dddm::pert_pep_small_temp" );
RealOptionKey const pert_pep_small_H( "dddm::pert_pep_small_H" );
RealOptionKey const pert_pep_small_L( "dddm::pert_pep_small_L" );
RealOptionKey const pert_pep_small_E( "dddm::pert_pep_small_E" );
RealOptionKey const pert_pep_shear_temp( "dddm::pert_pep_shear_temp" );
RealOptionKey const pert_pep_shear_H( "dddm::pert_pep_shear_H" );
RealOptionKey const pert_pep_shear_L( "dddm::pert_pep_shear_L" );
RealOptionKey const pert_pep_shear_E( "dddm::pert_pep_shear_E" );
RealOptionKey const pert_ter_small_temp( "dddm::pert_ter_small_temp" );
RealOptionKey const pert_ter_small_H( "dddm::pert_ter_small_H" );
RealOptionKey const pert_ter_small_L( "dddm::pert_ter_small_L" );
RealOptionKey const pert_ter_small_E( "dddm::pert_ter_small_E" );
RealOptionKey const pert_ter_shear_temp( "dddm::pert_ter_shear_temp" );
RealOptionKey const pert_ter_shear_H( "dddm::pert_ter_shear_H" );
RealOptionKey const pert_ter_shear_L( "dddm::pert_ter_shear_L" );
RealOptionKey const pert_ter_shear_E( "dddm::pert_ter_shear_E" );
IntegerOptionKey const pert_pep_num_rep( "dddm::pert_pep_num_rep" );
IntegerOptionKey const pert_ter_num_rep( "dddm::pert_ter_num_rep" );
IntegerOptionKey const pert_num( "dddm::pert_num" );
IntegerOptionKey const inner_num( "dddm::inner_num" );
RealOptionKey const ia_ener( "dddm::ia_ener" );

// design options
RealOptionKey const desn_mc_temp( "dddm::desn_mc_temp" );
}

// mover wrapper class
class DougsDockDesignMinimizeMagicMover : public Mover {

public:
	// default ctor
	DougsDockDesignMinimizeMagicMover(): Mover("DougsDockDesignMinimizeMagicMover"){}

	// dtor
	virtual ~DougsDockDesignMinimizeMagicMover(){}

	// methods
	void setup_pert_foldtree( core::pose::Pose & pose );
	void setup_filter_stats();
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "DougsDockDesignMinimizeMagicMover"; }

};

typedef utility::pointer::owning_ptr< DougsDockDesignMinimizeMagicMover > DougsDockDesignMinimizeMagicMoverOP;
typedef utility::pointer::owning_ptr< DougsDockDesignMinimizeMagicMover const > DougsDockDesignMinimizeMagicMoverCOP;

int
main( int argc, char* argv[] )
{
	try {
	/*********************************************************************************************************************
	  ____                                         ____  _       	  __  __
	 / ___|___  _ __ ___  _	__ ___ 	 ___  _	__		/ ___||	|_ _   _ / _|/ _|
	| |   /	_ \| '_	` _ \| '_ ` _ \	/ _ \| '_ \		\___ \|	__| | |	| |_| |_
	| |__| (_) | | | | | | | | | | | (_) | | | |	 ___) |	|_| |_|	|  _|  _|
	 \____\___/|_| |_| |_|_| |_| |_|\___/|_| |_|	|____/ \__|\__,_|_| |_|

	**********************************************************************************************************************/

	// add application specific options to options system
	option.add( dddm::pert_mc_temp, "The temperature to use for the pertubation phase of the DDDM protocol. Defaults to 0.8." ).def( 0.8 );
	option.add( dddm::pert_dock_rot_mag, "The rotation magnitude for the ridged body pertubation in the pertubation phase of the DDDM protocol. Defaults to 0.8." ).def( 0.5 );
	option.add( dddm::pert_dock_trans_mag, "The translation magnitude for the ridged body pertubation in the pertubation phase of the DDDM protocol. Defaults to 0.8." ).def( 0.25 );
	option.add( dddm::pert_pep_small_temp, "" ).def( 0.8 );
	option.add( dddm::pert_pep_shear_temp, "" ).def( 0.8 );
	option.add( dddm::pert_ter_small_temp, "" ).def( 0.8 );
	option.add( dddm::pert_ter_shear_temp, "" ).def( 0.8 );
	option.add( dddm::pert_pep_small_H, "" ).def( 1.0 );
	option.add( dddm::pert_pep_small_L, "" ).def( 1.0 );
	option.add( dddm::pert_pep_small_E, "" ).def( 1.0 );
	option.add( dddm::pert_pep_shear_H, "" ).def( 1.0 );
	option.add( dddm::pert_pep_shear_L, "" ).def( 1.0 );
	option.add( dddm::pert_pep_shear_E, "" ).def( 1.0 );
	option.add( dddm::pert_ter_small_H, "" ).def( 1.0 );
	option.add( dddm::pert_ter_small_L, "" ).def( 1.0 );
	option.add( dddm::pert_ter_small_E, "" ).def( 1.0 );
	option.add( dddm::pert_ter_shear_H, "" ).def( 1.0 );
	option.add( dddm::pert_ter_shear_L, "" ).def( 1.0 );
	option.add( dddm::pert_ter_shear_E, "" ).def( 1.0 );
	option.add( dddm::pert_pep_num_rep, "Number of small and shear iterations for the peptide" ).def( 100 );
	option.add( dddm::pert_ter_num_rep, "Number of small and shear iterations for the terminus" ).def( 100 );
	option.add( dddm::pert_num, "Number of iterations of perturbation loop per design" ).def(100);
	option.add( dddm::inner_num, "Number of iterations of the inner loop" ).def(100);
	option.add( dddm::ia_ener, "Upper energy limit for final design/interface analysis checkpoint" ).def( 0.0 );

	option.add( dddm::desn_mc_temp, "The temperature to use for the design/minimization phase of the DDDM protocol. Defaults to 0.8." ).def( 0.8 );

	// init command line options
	devel::init(argc, argv);

	// create an instance of my mover
	DougsDockDesignMinimizeMagicMoverOP D3DM( new DougsDockDesignMinimizeMagicMover() );

	// setup_filters
	D3DM->setup_filter_stats();

	// create job distributor
	protocols::jd2::JobDistributor::get_instance()->go( D3DM );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void
DougsDockDesignMinimizeMagicMover::apply(
   core::pose::Pose & pose
)
{

	// create score function
	TR << "Creating ScoreFunction..." << std::endl;
	core::scoring::ScoreFunctionOP score_fxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::MM_STD_WTS ) );


	// get a fold tree hopfully suitable for docking and loop modeling
	TR << "Creating Foldtree" << std::endl;
	setup_pert_foldtree( pose );

	// output starting pdb file
	//TR << "Output starting pdb..." << std::endl;
	//pose.dump_scored_pdb( "starting.pdb", *score_fxn );

	/*********************************************************************************************************************
	 ____  	       	_      	       _       	   _   _                ____  _
	|  _ \ ___ _ __| |_ _  	_ _ __|	|__   __ _| |_(_) ___  _ __		 |  _ \| |__   __ _ ___  ___
	| |_) /	_ \ '__| __| | | | '__|	'_ \ / _` | __|	|/ _ \|	'_ \	 | |_) |	_ \ / _` / __|/ _ \
	|  __/ 	__/ |  | |_| |_| | |  |	|_) | (_| | |_|	| (_) |	| | |	 |  __/| | | | (_| \__ \  __/
	|_|   \___|_|  	\__|\__,_|_|  |_.__/ \__,_|\__|_|\___/|_| |_|	 |_|   |_| |_|\__,_|___/\___|

	**********************************************************************************************************************/

	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, option[ dddm::pert_mc_temp ].value() ) );

	/*********************************************************
	  ___	       	 _   _ 	       	   ___ 	    _
	 |   \ ___  __| |_(_)_ _  __ _  / __|	___| |_	_  _ _ __
	 | |)	/ _ \/ _| / / |	' \/ _`	| \__ \/ -_)  _| || | '_ \
	 |___/\___/\__|_\_\_|_||_\__,	| |___/\___|\__|\_,_| .__/
	     	       	       	   |___/       	       	    |_|
	**********************************************************/

	// create a rigid body mover to move the peptide around in the pocket
	TR << "Setting up docking movers..." << std::endl;
	rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, option[ dddm::pert_dock_rot_mag].value(),  option[ dddm::pert_dock_trans_mag].value()) );
	//pert_dock_rbpm->rot_magnitude( option[ dddm::pert_dock_rot_mag].value() );
	//pert_dock_rbpm->trans_magnitude( option[ dddm::pert_dock_trans_mag].value() );


	/*********************************************************
	  ___	       	 _   _ 	  _    	  ___  	   _
	 | _ \___ _ __| |_(_)__| |___	 / __| ___| |_ _  _ _ __
	 |  _/ -_) '_	\  _| /	_` / -_) \__ \/	-_)  _|	|| | '_	\
	 |_| \___| .__/\__|_\__,_\___| |___/\___|\__|\_,_| .__/
	     	   |_| 	       	       	       	       	   |_|
	**********************************************************/
	TR << "Setting up peptide movers..." << std::endl;

	// get peptide start and end positions
	Size pep1_start( pose.conformation().chain_begin( 2 ) );
	Size pep1_end( pose.conformation().chain_end( 2 ) );

	// create movemap for peptide
	kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
	pert_pep_mm->set_bb_true_range(pep1_start, pep1_end);

	// create small and shear movers
	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, option[ dddm::pert_pep_small_temp ].value(), 1 ) );
	pert_pep_small->set_preserve_detailed_balance( true ); // skips the rama check in small move
	pert_pep_small->angle_max( 'H', option[ dddm::pert_pep_small_H ].value() );
	pert_pep_small->angle_max( 'L', option[ dddm::pert_pep_small_L ].value() );
	pert_pep_small->angle_max( 'E', option[ dddm::pert_pep_small_E ].value() );

	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( pert_pep_mm, option[ dddm::pert_pep_shear_temp ].value(), 1 ) );
	pert_pep_shear->set_preserve_detailed_balance( true ); // skips the rama check in shear move
	pert_pep_shear->angle_max( 'H', option[ dddm::pert_pep_shear_H ].value() );
	pert_pep_shear->angle_max( 'L', option[ dddm::pert_pep_shear_L ].value() );
	pert_pep_shear->angle_max( 'E', option[ dddm::pert_pep_shear_E ].value() );

	// create random mover
	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	pert_pep_random->add_mover( pert_pep_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, option[ dddm::pert_pep_num_rep ].value() ) );

	/*********************************************************
	  _____      	       _      _	  ___  	   _
	 |_  	_|__ _ _ _ __ (_)_ _ (_) / __| ___| |_ _  _ _ __
	   | |/ -_) '_| '  \|	| ' \| | \__ \/	-_)  _|	|| | '_	\
	   |_|\___|_|	|_|_|_|_|_||_|_| |___/\___|\__|\_,_| .__/
	     	       	       	       	       	       	   |_|
	**********************************************************/

	// get termini start and end positions
	Size ter1_start( pose.conformation().chain_begin( 1 ) );
	Size ter1_end( pose.conformation().chain_begin( 1 ) + 45 );

	// create movemap for termini
	kinematics::MoveMapOP pert_ter_mm( new kinematics::MoveMap() );
	pert_ter_mm->set_bb_true_range(ter1_start, ter1_end);

	// create small and shear movers
	simple_moves::SmallMoverOP pert_ter_small( new simple_moves::SmallMover( pert_ter_mm, option[ dddm::pert_ter_small_temp ].value(), 1 ) );
	pert_ter_small->set_preserve_detailed_balance( true ); // skips rama check in small move
	pert_ter_small->angle_max( 'H', option[ dddm::pert_ter_small_H ].value() );
	pert_ter_small->angle_max( 'L', option[ dddm::pert_ter_small_L ].value() );
	pert_ter_small->angle_max( 'E', option[ dddm::pert_ter_small_E ].value() );

 	simple_moves::ShearMoverOP pert_ter_shear( new simple_moves::ShearMover( pert_ter_mm,  option[ dddm::pert_ter_shear_temp ].value(), 1 ) );
	pert_ter_shear->set_preserve_detailed_balance( true ); // skips rama check in shear move
	pert_ter_shear->angle_max( 'H', option[ dddm::pert_ter_shear_H ].value() );
	pert_ter_shear->angle_max( 'L', option[ dddm::pert_ter_shear_L ].value() );
	pert_ter_shear->angle_max( 'E', option[ dddm::pert_ter_shear_E ].value() );

	// create random mover
	moves::RandomMoverOP pert_ter_random( new moves::RandomMover() );
	pert_ter_random->add_mover( pert_ter_small, 1 );
	pert_ter_random->add_mover( pert_ter_shear, 1 );

	// create repeat mover
	moves::RepeatMoverOP pert_ter_repeat( new moves::RepeatMover( pert_ter_random, option[ dddm::pert_ter_num_rep ].value() ) );


	/******************************************************************************
	  ___	    _  	       	       	  _____	   _   	  _    	 ___   	  _
	 | _ \___| |_	__ _ _ __  ___ _ |_   _| _(_)__	_| |___	/ __| ___| |_ _	 _ _ __
	 |   / _ \  _/ _` | '	 \/ -_)	'_|| ||	'_| / _` | (_-<	\__ \/ -_)  _| || | '_ \
	 |_|_\___/\__\__,_|_|_|_\___|_|  |_||_| |_\__,_|_/__/	|___/\___|\__|\_,_| .__/
	     	       	       	       	       	       	       	       	       	  |_|
	*******************************************************************************/
	TR << "Setting up RT movers..." << std::endl;

	// create a task factory and tack operations
	TaskFactoryOP pert_tf(new TaskFactory());

	operation::InitializeFromCommandlineOP pert_ifcl( new operation::InitializeFromCommandline() );
	pert_tf->push_back( pert_ifcl );

	operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	// create a rotamer trials mover
	protocols::simple_moves::RotamerTrialsMoverOP pert_rt(new protocols::simple_moves::EnergyCutRotamerTrialsMover( score_fxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	/*********************************************************
	   ___       	       	       	    ___	     _
	  / __|___ _ __  _ __	 ___ _ _   / __| ___| |_ _  _ _	__
	 | (__/ _ \ '	 \| '  \/ _ \ '	\  \__ \/ -_)  _| || | '_ \
	  \___\___/_|_|_|_|_|_\___/_||_| |___/\___|\__|\_,_| .__/
	     	       	       	       	       	       	     |_|
	**********************************************************/

	// create a random mover to hold the docking, termini, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 1 );
	//pert_random->add_mover( pert_ter_repeat, 1 );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	/*********************************************************************************************************************
	 ____  	       	 _		             __  __ _	         ____  _
	|  _ \ 	___  ___(_) __ _ _ __	    |  \/	 (_)_ __    |  _ \| |__	  __ _ ___  ___
	| | | |/ _ \/ __| |/ _`	| '_ \	  | |\/| | | '_	\   | |_) | '_ \ / _` /	__|/ _ \
	| |_| |	 __/\__	\ | (_|	| | | |	  | |  | | | | | |  |  __/| | |	| (_| \__ \  __/
	|____/ \___||___/_|\__,	|_| |_|	  |_|  |_|_|_| |_|  |_|	  |_| |_|\__,_|___/\___|
	       	       	   |___/
	**********************************************************************************************************************/
	/*********************************************************
	  ___	       	_      	      ___      _
	 |   \ ___ __(_)__ _ _ _   / __| ___|	|_ _  _	_ __
	 | |)	/ -_|_-< / _` |	' \  \__ \/ -_)	 _| || | '_ \
	 |___/\___/__/_\__, |_||_| |___/\___|\__|\_,_| .__/
	     	       	 |___/ 	       	       	       |_|
	**********************************************************/
	TR << "Setting up design movers..." << std::endl;

	// create a tack factory and task operations
	TaskFactoryOP desn_tf( new TaskFactory() );

	operation::InitializeFromCommandlineOP desn_ifc( new operation::InitializeFromCommandline() );
	desn_tf->push_back( desn_ifc );

	operation::ReadResfileOP desn_rrop( new operation::ReadResfile() );
	desn_rrop->default_filename();
	desn_tf->push_back( desn_rrop );

	// create a pack rotamers mover
	protocols::simple_moves::PackRotamersMoverOP desn_pr( new protocols::simple_moves::PackRotamersMover() );
	desn_pr->task_factory( desn_tf );
	desn_pr->score_function( score_fxn );
	desn_pr->nloop( 1 );


	/*********************************************************
	  __ 	__ _   	  _    	  _    	     ___      _
	 |  \/  (_)_ _ (_)_ __ (_)______  / __| ___| |_ _  _ _ __
	 | |\/| | | '	\| | ' 	\| |_ /	-_) \__	\/ -_) 	_| || |	'_ \
	 |_| 	|_|_|_||_|_|_|_|_|_/__\___| |___/\___|\__|\_,_|	.__/
	     	       	       	       	       	       	      |_|
	**********************************************************/
	TR << "Setting up minimization movers..." << std::endl;

	// create move map for minimization
	kinematics::MoveMapOP desn_mm( new kinematics::MoveMap() );

	// make jump minimizeable
	desn_mm->set_jump( 1, true );

	// make all the residues on the peptide we are moving minimizable
	for( Size i = pep1_start; i <= pep1_end; ++i ) {
		desn_mm->set_bb( i, true );
	}

	// create minimization mover
	protocols::simple_moves::MinMoverOP desn_min( new protocols::simple_moves::MinMover( desn_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01,	true ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	protocols::simple_moves::TaskAwareMinMoverOP desn_ta_min = new protocols::simple_moves::TaskAwareMinMover( desn_min, desn_tf );

	/*********************************************************
	   ___       	       	       	    ___	     _
	  / __|___ _ __  _ __	 ___ _ _   / __| ___| |_ _  _ _	__
	 | (__/ _ \ '	 \| '  \/ _ \ '	\  \__ \/ -_)  _| || | '_ \
	  \___\___/_|_|_|_|_|_\___/_||_| |___/\___|\__|\_,_| .__/
	     	       	       	       	       	       	     |_|
	**********************************************************/

	// create a sequence mover to hold pack rotamers and minimization movers
	moves::SequenceMoverOP desn_sequence( new moves::SequenceMover() );
	desn_sequence->add_mover( desn_pr );
	desn_sequence->add_mover( desn_ta_min );

	/*********************************************************************************************************************
	 __  __	      _	       	_
	|  \/  | __ _(_)_ __   | |    ___   ___	 _ __
	| |\/| |/ _` | | '_ \  | |   / _ \ / _ \| '_ \
	| |  | | (_| | | | | | | |__| (_) | (_)	| |_) |
	|_|  |_|\__,_|_|_| |_| |_____\___/ \___/| .__/
	       	       	       	       	       	|_|
	**********************************************************************************************************************/
	TR << "Starting main loop..." << std::endl;

	TR << *(desn_tf->create_task_and_apply_taskoperations( pose )) << std::endl;

	// get instance of job for output
	protocols::jd2::JobOP job_me( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// inner loop
	for ( Size k = 1; k <= Size( option[ dddm::inner_num ].value() ); ++k ) {

		// pert loop
		pert_mc->reset(pose);
		for( Size j = 1; j <= Size( option[ dddm::pert_num ].value() ); ++j )	{
			pert_trial->apply( pose );
			job_me->add_string_real_pair( "ENERGY_PERT", (*score_fxn)(pose) );
		}
		pert_mc->recover_low( pose );

		if ( k % 10 == 0 ) {
			// get packer task from task factory
			PackerTaskOP final_desn_pt( *(desn_tf->create_task_and_apply_taskoperations( pose )) );

			// add extra chi and extra chi cut off to pt
			for ( Size i = 1; i <= pose.total_residue(); ++i ) {
				final_desn_pt->nonconst_residue_task( i ).or_ex1( true );
				final_desn_pt->nonconst_residue_task( i ).or_ex2( true );
				final_desn_pt->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
			}

			// create a pack rotamers mover for the final design
			protocols::simple_moves::PackRotamersMoverOP final_desn_pr( new protocols::simple_moves::PackRotamersMover(score_fxn, final_desn_pt, 10 ) );

			// design with final pr mover
			final_desn_pr->apply( pose );

			// final min (okay to use ta min here)
			desn_ta_min->apply( pose );
		}
		else {
			// design
			desn_sequence->apply( pose );
			job_me->add_string_real_pair( "ENERGY_DESN", (*score_fxn)(pose) );
		}
	}

	TR << "Ending main loop..." << std::endl;

	TR << "Checking pose energy..." << std::endl;

		// create  MetricValues
	basic::MetricValue< core::Real > mv_sasa_complex;
	basic::MetricValue< core::Real > mv_sasa_seperated;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_complex;
	basic::MetricValue< utility::vector1< core::Size > > mv_unsat_res_seperated;
	basic::MetricValue< core::Real > mv_pack_complex;
	basic::MetricValue< core::Real > mv_pack_seperated;
	core::Real energy_complex;
	core::Real energy_seperated;
	core::Real hbond_ener_sum_complex;
	core::Real hbond_ener_sum_seperated;

	// calc energy
	energy_complex = (*score_fxn)(pose);

	// if score is greater than 0 this pose is no good and we should do a fail/retry
	if ( energy_complex > option[ dddm::ia_ener ].value() ) {
		TR << "Energy, " << energy_complex << " greater than cutoff, setting status to FAIL_RETRY..." << std::endl;
		set_last_move_status(protocols::moves::FAIL_RETRY);
		return;
	}
	// if pose is less than cutoff do a large final design and run filters
	else {

		TR << "Energy less than cutoff, doing final design and running filters..." << std::endl;

		// get packer task from task factory
		PackerTaskOP final_desn_pt( *(desn_tf->create_task_and_apply_taskoperations( pose )) );

		// add extra chi and extra chi cut off to pt
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			final_desn_pt->nonconst_residue_task( i ).or_ex1( true );
			final_desn_pt->nonconst_residue_task( i ).or_ex2( true );
			final_desn_pt->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
		}

		// create a pack rotamers mover for the final design
		protocols::simple_moves::PackRotamersMoverOP final_desn_pr( new protocols::simple_moves::PackRotamersMover(score_fxn, final_desn_pt, 10 ) );
		//final_desn_pr->packer_task( final_desn_pt );
		//final_desn_pr->score_function( score_fxn );
		//final_desn_pr->nloop( 10 );

		// design with final pr mover
		final_desn_pr->apply( pose );

		// final min (okay to use ta min here)
		desn_ta_min->apply( pose );

		// make copy of pose to calc stats
		Pose stats_pose( pose );

		// complex stats
		energy_complex = (*score_fxn)(stats_pose);
		stats_pose.metric("sasa","total_sasa",mv_sasa_complex);
		stats_pose.metric("unsat", "residue_bur_unsat_polars", mv_unsat_res_complex);
		utility::vector1< core::Size > const unsat_res_complex(mv_unsat_res_complex.value());
		stats_pose.metric( "pack", "total_packstat", mv_pack_complex );
		scoring::EnergyMap complex_emap( stats_pose.energies().total_energies() );
		hbond_ener_sum_complex = complex_emap[ hbond_sr_bb ] + complex_emap[ hbond_lr_bb ] + complex_emap[ hbond_bb_sc ] + complex_emap[ hbond_sc ];

		// seperate designed chain from other chains
		protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( pose, 1 ) ); // HARDCODED JUMP NUMBER
		translate->step_size( 1000.0 );
		translate->apply( stats_pose );

		// seperate stats
		energy_seperated = (*score_fxn)(stats_pose);
		stats_pose.metric("sasa","total_sasa",mv_sasa_seperated);
		stats_pose.metric("unsat", "residue_bur_unsat_polars", mv_unsat_res_seperated);
		utility::vector1< core::Size > const unsat_res_seperated(mv_unsat_res_seperated.value());
		stats_pose.metric( "pack", "total_packstat", mv_pack_seperated );
		scoring::EnergyMap seperated_emap( stats_pose.energies().total_energies() );
		hbond_ener_sum_seperated = seperated_emap[ hbond_sr_bb ] + seperated_emap[ hbond_lr_bb ] + seperated_emap[ hbond_bb_sc ] + seperated_emap[ hbond_sc ];

		// add values to job so that they will be output in the pdb
		job_me->add_string_real_pair( "ENERGY_COMPLEX:\t\t", energy_complex );
		job_me->add_string_real_pair( "ENERGY_SEPERATE:\t\t", energy_seperated );
		job_me->add_string_real_pair( "ENERGY_DIFF:\t\t", energy_complex - energy_seperated );

		job_me->add_string_real_pair( "SASA_COMPLEX:\t\t", mv_sasa_complex.value() );
		job_me->add_string_real_pair( "SASA_SEPERATE:\t\t", mv_sasa_seperated.value() );
		job_me->add_string_real_pair( "SASA_DIFF:\t\t", mv_sasa_complex.value() - mv_sasa_seperated.value() );

		job_me->add_string_real_pair( "HB_ENER_COMPLEX:\t\t", hbond_ener_sum_complex );
		job_me->add_string_real_pair( "HB_ENER_SEPERATE:\t\t", hbond_ener_sum_seperated );
		job_me->add_string_real_pair( "HB_ENER_DIFF:\t\t", hbond_ener_sum_complex - hbond_ener_sum_seperated );

		job_me->add_string_real_pair( "PACK_COMPLEX:\t\t", mv_pack_complex.value() );
		job_me->add_string_real_pair( "PACK_SEPERATE:\t\t", mv_pack_seperated.value() );
		job_me->add_string_real_pair( "PACK_DIFF:\t\t", mv_pack_complex.value() - mv_pack_seperated.value() );

	}

}

// this function setup the fold tree properly for docking
void
DougsDockDesignMinimizeMagicMover::setup_pert_foldtree(
  core::pose::Pose & pose
)
{
	using namespace kinematics;
	using namespace protocols::loops;
	using namespace std;

	// get current fold tree
	FoldTree f( pose.fold_tree() );
	f.clear();

	/*********************************************************
	  ___	       	 _     	  _
	 |   \ ___  __| |__  _ | |_  _ _ __  _ __ ___
	 | |)	/ _ \/ _| / / |	|| | ||	| '  \|	'_ (_-<
	 |___/\___/\__|_\_\  \__/ \_,_|_|_|_|	.__/__/
	     	       	       	       	      |_|
	**********************************************************/

	// get the start and end for both chains
	Size pro_start( pose.conformation().chain_begin( 1 ) );
	Size pro_end( pose.conformation().chain_end( 1 ) );
	Size pep_start( pose.conformation().chain_begin( 2 ) );
	Size pep_end( pose.conformation().chain_end( 2 ) );

	// get jump positions based on the center of mass of the chains
	Size dock_jump_pos_pro( core::pose::residue_center_of_mass( pose, pro_start, pro_end ) );
	Size dock_jump_pos_pep( core::pose::residue_center_of_mass( pose, pep_start, pep_end ) );

	// build fold tree
	Size jump_index( f.num_jump() + 1 );
	f.add_edge( pro_start, dock_jump_pos_pro, -1 );
	f.add_edge( dock_jump_pos_pro, pro_end, -1 );
	f.add_edge( pep_start, dock_jump_pos_pep, -1);
	f.add_edge( dock_jump_pos_pep, pep_end, -1 );
	f.add_edge( dock_jump_pos_pro, dock_jump_pos_pep, jump_index );

	// set pose foldtree to foldtree we just created
	f.reorder(1);
	f.check_fold_tree();
	assert( f.check_fold_tree() );

	std::cout << "FOLDTREE: " << f << std::endl;

	pose.fold_tree( f );

	// DEBUG
	//core::kinematics::dump_pose_kinemage( "pose.kin", pose );

}

void
DougsDockDesignMinimizeMagicMover::setup_filter_stats()
{
	/*********************************************************************************************************************
  _____	_ _ _  	       	      __  ____ 	_      	 _           ____	      _
 |  ___(_) | |_	___ _ __     / / / ___|| |_ __ _| |_ ___	  / ___|  ___| |_ _   _	_ __
 | |_  | | | __/ _ \ '__|   / /	 \___ \| __/ _`	| __/ __|	  \___ \ / _ \ __| | | | '_ \
 |  _| | | | ||	 __/ | 	   / / 	  ___) | || (_|	| |_\__	\	   ___)	|  __/ |_| |_| | |_) |
 |_|   |_|_|\__\___|_| 	  /_/  	 |____/	\__\__,_|\__|___/	  |____/ \___|\__|\__,_| .__/
																													       	       	       |_|
	*********************************************************************************************************************/

	// create and register sasa calculator
	pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

	// create and register hb calculator
	pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );

	// create and register unsat calculator
	pose::metrics::PoseMetricCalculatorOP unsat_calculator( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds") ) ;
	pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );

	// create and register packstat calculator
	pose::metrics::PoseMetricCalculatorOP pack_calcculator( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( "pack", pack_calcculator );

}
