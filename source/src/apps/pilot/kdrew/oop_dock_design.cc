// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//kdrew: based on code from mini/src/apps/public/scenarios/doug_dock_design_min_mod2_cal_cal.cc
//   and https://svn.rosettacommons.org/source/branches/releases/rosetta-3.1/manual/advanced/example_protocol.cc

//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers

#include <protocols/ncbb/oop/OopDockDesignProtocol.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/ncbb/oop/OopRandomPuckMover.hh>
#include <protocols/ncbb/oop/OopRandomSmallMover.hh>
#include <protocols/ncbb/oop/OopPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::ncbb::oop;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::oop;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// Here is a sketch of the basic flow of the program...
//
// Pertubation Phase
//   +-Monte Carlo Mover---------------------------------------+
//   | +-Random Mover-( 1 / 2 / 1 / 1 )--------------------+ | |
//   | | +-Docking Mover-----------------------------------+ | |
//   | | | small rigid body movements between the peptide  | | |
//   | | | and protein for conformational diversity        | | |
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
//   | +-Pack Rotamers Mover---------------------------------+ |
//   | | repack and design rotamers to explore sequence     | |
//   | | space                                              | |
//   | +-----------------------------------------------------+ |
//   | +-Minimization Mover----------------------------------+ |
//   | | energy minimize the current conformation            | |
//   | +-----------------------------------------------------+ |


// tracer - used to replace cout
static basic::Tracer TR( "ODDM" );

// application specific options
namespace oddm {
// pert options
RealOptionKey const mc_temp( "oddm::mc_temp" );
RealOptionKey const pert_mc_temp( "oddm::pert_mc_temp" );
RealOptionKey const pert_dock_rot_mag( "oddm::pert_dock_rot_mag" );
RealOptionKey const pert_dock_trans_mag( "oddm::pert_dock_trans_mag" );
RealOptionKey const pert_pep_small_temp( "oddm::pert_pep_small_temp" );
RealOptionKey const pert_pep_small_H( "oddm::pert_pep_small_H" );
RealOptionKey const pert_pep_small_L( "oddm::pert_pep_small_L" );
RealOptionKey const pert_pep_small_E( "oddm::pert_pep_small_E" );
RealOptionKey const pert_pep_shear_temp( "oddm::pert_pep_shear_temp" );
RealOptionKey const pert_pep_shear_H( "oddm::pert_pep_shear_H" );
RealOptionKey const pert_pep_shear_L( "oddm::pert_pep_shear_L" );
RealOptionKey const pert_pep_shear_E( "oddm::pert_pep_shear_E" );

IntegerOptionKey const pert_pep_num_rep( "oddm::pert_pep_num_rep" );
IntegerOptionKey const pert_num( "oddm::pert_num" );
IntegerOptionKey const dock_design_loop_num( "oddm::dock_design_loop_num" );

BooleanOptionKey const no_design( "oddm::no_design" );
BooleanOptionKey const final_design_min( "oddm::final_design_min" );
BooleanOptionKey const use_soft_rep( "oddm::use_soft_rep" );
BooleanOptionKey const mc_initial_pose( "oddm::mc_initial_pose" );
BooleanOptionKey const oop_design_first( "oddm::oop_design_first" );

BooleanOptionKey const pymol( "oddm::pymol" );
BooleanOptionKey const keep_history( "oddm::keep_history" );

// design options
RealOptionKey const desn_mc_temp( "oddm::desn_mc_temp" );

//IntegerVectorOptionKey const oop_positions( "oddm::oop_positions" );

}

/*
* class OopDockDesignMinimizeMover : public Mover {

public:

//default ctor
OopDockDesignMinimizeMover(): Mover("OopDockDesignMinimizeMover"){}

//default dtor
virtual ~OopDockDesignMinimizeMover(){}

//methods
void setup_pert_foldtree( core::pose::Pose & pose);
void setup_filter_stats();
virtual void apply( core::pose::Pose & pose );
virtual std::string get_name() const { return "OopDockDesignMinimizeMover"; }

};

typedef utility::pointer::owning_ptr< OopDockDesignMinimizeMover > OopDockDesignMinimizeMoverOP;
typedef utility::pointer::owning_ptr< OopDockDesignMinimizeMover const > OopDockDesignMinimizeMoverCOP;
*/

int
main( int argc, char* argv[] )
{
	try {
		/*********************************************************************************************************************
		Common Setup
		**********************************************************************************************************************/

		// add application specific options to options system
		// There are far more options here than you will realistically need for a program of this complexity - but this gives you an idea of how to fine-grain option-control everything
		option.add( oddm::mc_temp, "The temperature to use for the outer loop of the ODDM protocol. Defaults to 1.0." ).def( 1.0 );
		option.add( oddm::pert_mc_temp, "The temperature to use for the pertubation phase of the ODDM protocol. Defaults to 0.8." ).def( 0.8 );
		option.add( oddm::pert_dock_rot_mag, "The rotation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 1.0." ).def( 1 );
		option.add( oddm::pert_dock_trans_mag, "The translation magnitude for the ridged body pertubation in the pertubation phase of the ODDM protocol. Defaults to 0.5." ).def( 0.5 );
		option.add( oddm::pert_pep_small_temp, "" ).def( 0.8 );
		option.add( oddm::pert_pep_shear_temp, "" ).def( 0.8 );

		option.add( oddm::pert_pep_small_H, "" ).def( 2.0 );
		option.add( oddm::pert_pep_small_L, "" ).def( 2.0 );
		option.add( oddm::pert_pep_small_E, "" ).def( 2.0 );
		option.add( oddm::pert_pep_shear_H, "" ).def( 2.0 );
		option.add( oddm::pert_pep_shear_L, "" ).def( 2.0 );
		option.add( oddm::pert_pep_shear_E, "" ).def( 2.0 );

		option.add( oddm::pert_pep_num_rep, "Number of small and shear iterations for the peptide" ).def( 100 );
		option.add( oddm::pert_num, "Number of iterations of perturbation loop per design" ).def(10);
		option.add( oddm::dock_design_loop_num, "Number of iterations of pertubation and design" ).def(10);

		option.add( oddm::no_design, "Only repack, do not design. Default false" ).def(false);
		option.add( oddm::final_design_min, "Do a final repack/design and minimization. Default true" ).def(true);
		option.add( oddm::use_soft_rep, "Use soft repulsion for pertubation and initial design. Default false" ).def(false);
		option.add( oddm::mc_initial_pose, "Allow initial pose to be considered as lowest energy pose. Default false" ).def(false);
		option.add( oddm::oop_design_first, "Design before pertubation (want when initial struct is aligned to hotspot)  Default false" ).def(false);

		option.add( oddm::pymol, "Set up pymol mover. Default false" ).def(false);
		option.add( oddm::keep_history, "Keep history in pymol. Requires oddm::pymol set to true. Default false" ).def(false);

		option.add( oddm::desn_mc_temp, "The temperature to use for the design/minimization phase of the ODDM protocol. Defaults to 0.8." ).def( 0.8 );

		//utility::vector1< core::Size > empty_vector(0);
		//option.add( oddm::oop_positions, "The positions of the first residues of oop rings" ).def( empty_vector );

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		scoring::ScoreFunctionOP score_fxn = get_score_function();
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

		//create mover instance
		OopDockDesignProtocolOP ODDM_mover( new OopDockDesignProtocol(
			score_fxn,
			option[ oddm::mc_temp].value(),
			option[ oddm::pert_mc_temp].value(),
			option[ oddm::pert_dock_rot_mag].value(),
			option[ oddm::pert_dock_trans_mag].value(),
			option[ oddm::pert_pep_small_temp].value(),
			option[ oddm::pert_pep_small_H].value(),
			option[ oddm::pert_pep_small_L].value(),
			option[ oddm::pert_pep_small_E].value(),
			option[ oddm::pert_pep_shear_temp].value(),
			option[ oddm::pert_pep_shear_H].value(),
			option[ oddm::pert_pep_shear_L].value(),
			option[ oddm::pert_pep_shear_E].value(),
			option[ oddm::pert_pep_num_rep].value(),
			option[ oddm::pert_num].value(),
			option[ oddm::dock_design_loop_num].value(),
			option[ oddm::no_design].value(),
			option[ oddm::final_design_min].value(),
			option[ oddm::use_soft_rep].value(),
			option[ oddm::mc_initial_pose].value(),
			option[ oddm::oop_design_first].value(),
			option[ oddm::pymol].value(),
			option[ oddm::keep_history].value()
			) );

		protocols::ncbb::setup_filter_stats();

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( ODDM_mover );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main
