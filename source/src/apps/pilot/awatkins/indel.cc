// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//awatkins: based heavily on code from kdrew/oop_dock_design.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMoverCreator.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>


#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <core/fragment/FragSet.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/make_loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loop_build/LoopBuildMover.hh>


#include <protocols/comparative_modeling/LoopRelaxMover.hh>

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
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

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
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// Here is a sketch of the basic flow of the program...
//
// Pertubation Phase
//   +-Monte Carlo Mover---------------------------------------+
//   | +-Random Mover-( 1 / 2 / 1 / 1 )--------------------+ | |
//   | | +-Docking Mover-----------------------------------+ | |
//   | | | small rigid body movements between the peptide	 | | |
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
//   | | repack and design rotamers to explore sequence 	   | |
//   | | space  	                                           | |
//   | +-----------------------------------------------------+ |
//   | +-Minimization Mover----------------------------------+ |
//   | | energy minimize the current conformation            | |
//   | +-----------------------------------------------------+ |


// tracer - used to replace cout
static basic::Tracer TR("IndelOptimization`");

// application specific options
namespace indel {
	// pert options
	IntegerOptionKey const start_res( "indel::start_res" );
	IntegerOptionKey const end_res( "indel::end_res" );
	IntegerOptionKey const loop_length( "indel::loop_length" );
}

class IndelOptimizationMover : public Mover {

	public:

		//default ctor
		IndelOptimizationMover(): Mover("IndelOptimizationMover"){}

		//default dtor
		virtual ~IndelOptimizationMover(){}

		//methods
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const { return "IndelOptimizationMover"; }

};

typedef utility::pointer::shared_ptr< IndelOptimizationMover > IndelOptimizationMoverOP;
typedef utility::pointer::shared_ptr< IndelOptimizationMover const > IndelOptimizationMoverCOP;


int
main( int argc, char* argv[] )
{
try {
	/*********************************************************************************************************************
	Common Setup
	***************************( *******************************************************************************************/

	// add application specific options to options system
		
	option.add( indel::start_res, "The first residue to delete" ).def( 1 );
	option.add( indel::end_res, "The last residue to delete (will be set to start_res if not specified)" ).def( 1 );
	option.add( indel::loop_length, "The number of residues behind and in front of the deleted residue(s) to remodel" ).def( 4 );
		
	devel::init(argc, argv);

	//create mover instance
	IndelOptimizationMoverOP indel_mover( new IndelOptimizationMover() );

	//call job distributor
	protocols::jd2::JobDistributor::get_instance()->go( indel_mover );
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
}//main

void
IndelOptimizationMover::apply(
	core::pose::Pose & pose
)
{
	using namespace scoring;
	using namespace constraints;
	using namespace func;
	ScoreFunctionOP score_fxn = get_score_function();
	add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	add_fa_constraints_from_cmdline_to_pose(pose);

	Size start_res = option[ indel::start_res ].value();
	Size end_res;
	if ( option[ indel::end_res ].user() ) {
		end_res = option[ indel::end_res ].value();
	} else {
		end_res = start_res;
	}
	
	Size half_loop_length = option[ indel::loop_length ].value();
	
	pose.conformation().delete_residue_range_slow( start_res, end_res );
	pose.conformation().detect_disulfides();
	pose.conformation().declare_chemical_bond( start_res-1, "C", end_res, "N" );
	TR << pose.fold_tree();
	
	// Add constraint
	/*if ( resnum > 1 ) {
		TR << "adding AtomPairConstraint between residues " << resnum-1 << " and " << resnum << std::endl;
		
		HarmonicFuncOP harm_func( new HarmonicFunc( 1.33, .02 ) );
		AtomID aidC( pose.residue( resnum-1 ).atom_index("C"), resnum-1 );
		AtomID aidN( pose.residue( resnum ).atom_index("N"), resnum );
		ConstraintCOP atompair( ConstraintOP( new AtomPairConstraint( aidC, aidN, harm_func ) ) );
		pose.add_constraint( atompair );
		score_fxn->set_weight( atom_pair_constraint, 1 );
	}*/

	// Loop model the remaining segment; range size before and
	Size loop_start = ( start_res < 2 + half_loop_length ) ? 2 : start_res - half_loop_length;
	Size loop_end   = ( end_res + half_loop_length > pose.total_residue() ) ? pose.total_residue() : end_res + half_loop_length;
	
	TR << "loop is from " << loop_start << " to " << loop_end << std::endl;
	
	Size cutpoint = Size( ( loop_start + loop_end ) / 2 );
	if ( cutpoint == loop_start ) ++cutpoint;
	if ( cutpoint == loop_end   ) --cutpoint;
	//protocols::loops::Loop loop( loop_start, loop_end, cutpoint );
	//TR << "i.e. " << loop << std::endl;
	protocols::loops::LoopsOP loops( new protocols::loops::Loops );
	loops->push_back( loop_start, loop_end, cutpoint );//loop );
	TR << (*loops) << std::endl;
	
	// Add constraints to native except loop
	core::pose::PoseOP native_pose( new core::pose::Pose( pose ) );
	native_pose->conformation().delete_residue_range_slow(loop_start, loop_end);
	native_pose->conformation().detect_disulfides();
	
	protocols::relax::AtomCoordinateCstMover coord_cst;
	coord_cst.set_refstruct( native_pose );
	coord_cst.apply( pose );

	// Set up LoopRelaxMover
	std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
	std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
	std::string const refine       ( option[ OptionKeys::loops::refine ]() );
	std::string const relax        ( option[ OptionKeys::loops::relax ]() );
	
	// fragment initialization (not necessary in this case)
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( option[ OptionKeys::loops::frag_files ].user()) {
		protocols::loops::read_loop_fragments( frag_libs );
	}
	
	//setup of looprelax_mover
	protocols::comparative_modeling::LoopRelaxMover looprelax_mover;//( remodel, intermedrelax, refine, relax, loops);
	looprelax_mover.loops( loops );
	TR << (*looprelax_mover.get_loops() ) << std::endl;
	looprelax_mover.frag_libs( frag_libs );
	looprelax_mover.relax( relax );
	looprelax_mover.refine( refine );
	looprelax_mover.remodel( remodel );
	looprelax_mover.intermedrelax( intermedrelax );
	
	protocols::loop_build::LoopBuildMoverOP loopbuild_mover( new protocols::loop_build::LoopBuildMover(looprelax_mover) );
	loopbuild_mover->apply( pose );
	
	return;
}
