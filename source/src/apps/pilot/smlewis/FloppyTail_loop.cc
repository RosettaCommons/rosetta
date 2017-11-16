// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/scenarios/FloppyTail.cc
/// @brief This app is intended for modeling the binding of a long unstructured C-terminal tail to some other part of a protein.  It will work for internal regions too.
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>

#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/constraints/util.hh>

#include <protocols/moves/MonteCarlo.hh>

//loop-related machinery
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>


//movers
#include <protocols/simple_moves/BackboneMover.hh> //SmallMover
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> //typeset swapping
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/OutputMovers.hh> //pdbdumpmover

//calculators and neighbor detection machinery
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <protocols/moves/MoverStatistics.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.public.scenarios.FloppyTail" );

namespace basic { namespace options { namespace OptionKeys {
basic::options::IntegerOptionKey const flexible_start_resnum("flexible_start_resnum");
basic::options::IntegerOptionKey const flexible_stop_resnum("flexible_stop_resnum");
basic::options::StringOptionKey const flexible_chain("flexible_chain");
basic::options::RealOptionKey const shear_on("shear_on");
basic::options::RealOptionKey const short_tail_fraction("short_tail_fraction");
basic::options::RealOptionKey const short_tail_off("short_tail_off");
basic::options::BooleanOptionKey const pair_off("pair_off");
basic::options::BooleanOptionKey const publication("publication");
basic::options::BooleanOptionKey const C_root("C_root");
}}}//basic::options::OptionKeys

/// @details This function is specific to the original system for which this code was written - if you are not trying to duplicate the initial results you should remove it!
void create_extra_output( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP score_fxn ){

	//score main pose
	core::Real const full_score((*score_fxn)(pose));

	//constants for hacking up pose
	core::Size const E2_chain(pose.conformation().num_chains()); //E2 is last chain
	core::Size const E2_start(pose.conformation().chain_begin(E2_chain));

	//create and score E3_RING portion
	core::pose::Pose E3_RING(pose);
	E3_RING.conformation().delete_residue_range_slow(E2_start, E3_RING.size());
	core::Real const E3_RING_score((*score_fxn)(E3_RING));

	//create and score E2 portion
	core::pose::Pose E2 = pose.split_by_chain(E2_chain);
	core::Real const E2_score((*score_fxn)(E2));

	//print Job
	using protocols::jd2::JobDistributor;
	protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
	//debugging - check these pdbs
	//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, E2, "E2");
	//JobDistributor::get_instance()->job_outputter()->other_pose( job_me, E3_RING, "E3_RING");

	//print Binding_Energy
	core::Real const BE(full_score - E3_RING_score - E2_score);
	job_me->add_string("Binding_Energy = Complex - E2 alone - E3/RING alone");
	job_me->add_string_real_pair("Binding energy", BE);
	//print crosslink distance
	//magic numbers: crosslink experiment tested these residues
	bool const NEDD8(E2_chain == 4);
	core::Size const cl_base(pose.pdb_info()->pdb2pose('A', (NEDD8 ? 1676: 679)));   //K1676??
	core::Size const cl_tail(pose.pdb_info()->pdb2pose('C', 223));

	core::Real const cl_dist(pose.residue(cl_base).atom("CA").xyz().distance( pose.residue(cl_tail).atom("CA").xyz() ));
	job_me->add_string("Crosslink = distance between atom CA on residues (chain A, K679/K1676; chain C, 223)");
	job_me->add_string_real_pair("crosslink distance", cl_dist);

	//std::cout << full_score << " " << E3_RING_score << " " << E2_score << std::endl;

	return;
}

/// @brief FloppyTail mover
class FloppyTailMover : public protocols::moves::Mover {
public:
	FloppyTailMover() : start_(0), stop_(0), init_for_input_yet_(false), centroid_scorefunction_(NULL),
		fullatom_scorefunction_(NULL), task_factory_(NULL), movemap_(NULL), movemap_lesstail_(NULL),
		foldtree_(NULL), fragset3mer_(NULL), loops_()
	{
		//this should be per-input, not in the ctor, if multiple frags files
		if ( basic::options::option[ basic::options::OptionKeys::in::file::frag3].user() ) {
			fragset3mer_ = new core::fragment::ConstantLengthFragSet( 3 );
			fragset3mer_->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag3].value() );
		}

		bool const pair_off(basic::options::option[ basic::options::OptionKeys::pair_off ].value() );

		//set up centroid scorefunction
		using namespace core::scoring;
		centroid_scorefunction_ = new ScoreFunction;
		centroid_scorefunction_->set_weight( env,         1.0 );
		centroid_scorefunction_->set_weight( cbeta,       1.0 );
		centroid_scorefunction_->set_weight( vdw,         1.0 );
		centroid_scorefunction_->set_weight( pair, (pair_off ? 0.0 : 1.0) ); //no pair term experiment - not for general use
		centroid_scorefunction_->set_weight( cenpack,     1.0 );
		centroid_scorefunction_->set_weight( rama,        1.0 );
		centroid_scorefunction_->set_weight( hbond_lr_bb, 1.0 );
		centroid_scorefunction_->set_weight( hbond_sr_bb, 1.0 );
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *centroid_scorefunction_ ); //protected if(option) internally
		TR << "Using centroid scorefunction\n" << *centroid_scorefunction_;

		//set up fullatom scorefunction
		fullatom_scorefunction_ = get_score_function();
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ ); //protected if(option) internally
		if ( pair_off ) fullatom_scorefunction_->set_weight( fa_pair, 0.0 ); //not for general use
		TR << "Using fullatom scorefunction (TALARIS_2013), pair may be modified\n"
			<< *fullatom_scorefunction_;

	}

	/// @brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
	virtual
	void
	init_on_new_input(core::pose::Pose const & pose) {
		init_for_input_yet_ = true;

		//determine where the flexible tail is
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if ( !option[in::file::movemap].user() ) {
			//original code: read from command line options
			char const chain(option[flexible_chain].value()[0]); //just take the first one
			start_ = pose.pdb_info()->pdb2pose().find(chain, option[flexible_start_resnum].value());
			if ( option[flexible_stop_resnum].user() ) {
				stop_ = pose.pdb_info()->pdb2pose().find(chain, option[flexible_stop_resnum].value());
			} else {
				stop_ = pose.size();
			}
			TR << "Tail is from " << start_ << " to " << stop_ << std::endl;

			//setup MoveMap
			movemap_ = new core::kinematics::MoveMap;
			for ( core::Size i(start_); i<=stop_; ++i ) {
				movemap_->set_bb(i, true); //backbone mobile
				movemap_->set_chi(i, true); //chi of mobile residues
			}

			movemap_lesstail_ = new core::kinematics::MoveMap(*movemap_);
			if ( stop_ == pose.size() ) {
				core::Size const taillength = stop_ - start_;
				core::Size const substop = (start_ + core::Size(core::Real(taillength) * (1.0 - option[ short_tail_fraction ] ) ) );
				for ( core::Size i(start_); i <= substop; ++i ) {
					movemap_lesstail_->set_bb(i, false);
					movemap_lesstail_->set_chi(i, false);
				}
			}
		} else if ( //these options are incompatible with movemap use; the movemap determines flexible regions
				option[flexible_chain].user() ||
				option[flexible_stop_resnum].user() ||
				option[flexible_chain].user() ||
				option[short_tail_fraction].user() ||
				option[short_tail_off].user() ) {
			utility_exit_with_message("option in::file::movemap not compatible with options flexible_chain, flexible_stop_resnum, flexible_chain, short_tail_fraction, or short_tail off.  This is because a manually-defined movemap overrides these options.");
		} else {
			//handle user-defined movemap from file; reverse-convert into start_ and stop_
			movemap_ = new core::kinematics::MoveMap;
			movemap_->init_from_file(option[in::file::movemap].value());
			movemap_lesstail_ = new core::kinematics::MoveMap(*movemap_);

			//calculate effective start_ and stop_.  This may be less efficient than the MoveMap's iterators, but it is vastly simpler, less likely to be buggy, and not a performace concern.
			core::Size const nres(pose.size());
			for ( core::Size i(1); i<=nres; ++i ) {
				if ( movemap_->get_bb(i) ) {
					start_ = i;
					break;
				}
			}

			for ( core::Size i(nres); i>=1; --i ) {
				if ( movemap_->get_bb(i) ) {
					stop_ = i;
					break;
				}
			}

		}

		//error handle: if start_ and stop_ couldn't be set from movemap somehow, freak out
		if ( (start_ == 0) || (stop_ == 0) ) {
			std::ostringstream message;
			message << "invalid flexible region (start (" << start_ << ") or stop (" << stop_ << ") is undefined) - check your flags or movemap file";
			utility_exit_with_message( message.str());
		}

		//We want to linearize the fold_tree so that internal linkers with noncovalent attachments on both sides will work properly; it will have no effect on terminal tails
		foldtree_ = new core::kinematics::FoldTree(pose.fold_tree()); //store original fold tree; if we enter neither option below it stays valid
		if ( !((stop_ == pose.size()) || (start_ == 1)) ) {
			TR << "non-terminal or N-terminal (C-rooted) floppy section, using a linear fold tree to try to ensure downstream residues follow \nOld tree: " << pose.fold_tree();
			foldtree_ = new core::kinematics::FoldTree(core::kinematics::linearize_fold_tree(*foldtree_));
			TR << "new tree: " << *foldtree_ << std::endl;
		}
		if ( basic::options::option[ basic::options::OptionKeys::C_root].value() ) {
			foldtree_->reorder(pose.size());
			TR << "C-rooted tree: " << *foldtree_ << std::endl;
		}


		//setup of loops objects, as needed
		std::set< core::Size > loop_posns; //may be needed for neighbor detection
		if ( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user() ) {
			loops_ = protocols::loops::Loops( true ); //automatically reads that command line
			TR << loops_ << std::endl;
			//set up interface-plus-neighbors-positions operation
			for ( protocols::loops::Loops::const_iterator it(loops_.begin()), end(loops_.end()); it != loops_.end(); it++ ) {
				for ( core::Size j(it->start()), end(it->stop()); j <= end; ++j ) {
					loop_posns.insert(j); //track for packing
				}//for each residue in loop
			}//for each loop
		}//if loops exist


		//setup of TaskFactory
		//command line and resfile options
		using namespace core::pack::task;
		task_factory_ = new TaskFactory;
		task_factory_->push_back( new operation::InitializeFromCommandline );
		if ( option[ packing::resfile ].user() ) {
			task_factory_->push_back( new operation::ReadResfile );
		}

		//iterate through movemap, determining where regions of flexibility and inflexibility are
		//Each inflexible region surrounded by two flexible region constitutes a "group" for InterGroupNeighborsCalculator
		//inclusion of flexible regions into multiple groups ensures they will repack regardless of their neighbors
		//all group pairs are passed as viable to InterGroupNeighborsCalculator
		//flexible regions are included in two subsequent groups, ensuring that they'll pack at all times (as their own neighbors)
		bool previous_state(false); //assume we start with a non-flexible region; new regions are triggered on flexiblity
		utility::vector1< std::set < core::Size > > regions; //a set of regions to turn into groups for comparison
		std::set < core::Size > const empty; //easier to add empty sets to the vector than construct-then-add
		core::Size current_group(1);
		regions.push_back(empty);

		//iterate through all residues in the pose/movemap
		core::Size const nres(pose.size());
		for ( core::Size i(1); i<=nres; ++i ) {
			bool const this_state(movemap_->get_bb(i));
			if ( this_state != previous_state ) { //if we are changing regions
				if ( previous_state == false ) { //we are leaving an inflexible region
					regions.push_back(empty); //add a new group, and start putting residues there
					++current_group;
				}
				previous_state = this_state; //keep previous_state up to date
			}
			regions[current_group].insert(i); //add this residue to the current group
			if ( this_state ) regions[current_group-1].insert(i); //add this residue to the previous group, if flexible - this ensures flexible regions always can repack
		}

		//make all pairs of groups (without replacement
		//if you have 1, 2, 3, 4; make 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
		core::Size const num_regions(regions.size());
		utility::vector1< std::pair< std::set<core::Size>, std::set<core::Size> > > vector_of_pairs;
		for ( core::Size first_group(1); first_group < num_regions; ++first_group ) {
			for ( core::Size second_group(first_group+1); second_group <= num_regions; ++second_group ) {
				vector_of_pairs.push_back(std::make_pair(regions[first_group], regions[second_group]));
			}
		}

		//check contents of vector_of_pairs
		core::Size const num_pairs(vector_of_pairs.size());
		for ( core::Size i(1); i<=num_pairs; ++i ) {
			core::Size const
				onestart(*(vector_of_pairs[i].first.begin())),
				onestop(*(vector_of_pairs[i].first.rbegin())),
				twostart(*(vector_of_pairs[i].second.begin())),
				twostop(*(vector_of_pairs[i].second.rbegin()));

			TR << "IGNC will compare group " << onestart << "-" << onestop << " with " << twostart << "-" << twostop << std::endl;

			core::Size guess(onestart);
			for ( std::set<core::Size>::const_iterator iter(vector_of_pairs[i].first.begin()), end(vector_of_pairs[i].first.end()); iter != end; ++iter ) {
				if ( guess++ != *iter ) TR.Error << "non-contiguous set, debug me!" << std::endl;
				//TR << *iter << std::endl;
			}
			guess = twostart;
			for ( std::set<core::Size>::const_iterator iter(vector_of_pairs[i].second.begin()), end(vector_of_pairs[i].second.end()); iter != end; ++iter ) {
				if ( guess++ != *iter ) TR.Error << "non-contiguous set, debug me!" << std::endl;
				//TR << *iter << std::endl;
			}

		}

		//check if calculator exists; create if not
		std::string const calc("IGNC_FloppyTail");
		if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc) ) {
			core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc);
			TR << "removed a PoseMetricCalculator " << calc << ", hopefully this is due to multiple inputs to FloppyTail and not a name clash" << std::endl;
		}
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc, new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(vector_of_pairs) );

		//now that calculator exists, add the sucker to the TaskFactory via RestrictByCalculatorsOperation
		utility::vector1< std::pair< std::string, std::string > > calculators_used;
		std::pair< std::string, std::string > IGNC_cmd( calc, "neighbors" );
		calculators_used.push_back( IGNC_cmd );

		//also create calculator; etc for loops if they exist
		if ( !loops_.empty() && !loop_posns.empty() ) {
			TR << "also creating NeighborhoodByDistanceCalculator for loops" << std::endl;
			std::string const calc_loop("NBDC_FloppyTail");
			//check if it exists, remove if so
			if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc_loop) ) {
				core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc_loop);
				TR << "removed a PoseMetricCalculator " << calc_loop << ", hopefully this is due to multiple inputs to FloppyTail and not a name clash" << std::endl;
			}
			//create and register calculator
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc_loop, new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) );

			std::pair< std::string, std::string > NBDC_cmd( calc_loop, "neighbors" );
			calculators_used.push_back( NBDC_cmd );
		}

		task_factory_->push_back( new protocols::toolbox::task_operations::RestrictByCalculatorsOperation( calculators_used ) );


		//debugging: print PackerTask
		//TR << *(task_factory_->create_task_and_apply_taskoperations( pose )) << std::endl;

		return;
	}

	virtual ~FloppyTailMover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		if ( !init_for_input_yet_ ) init_on_new_input(pose);

		//apply fold tree (determined in init_on_new_input) and report
		pose.fold_tree(*foldtree_);
		TR << "foldtree, movemap: " << std::endl;
		core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, TR);
		if ( stop_ == pose.size() ) {
			TR << "foldtree, movemap for first part of refine: " << std::endl;
			core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_lesstail_, TR);
		}

		core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose ); //protected internally if no constraints

		//centroid
		clock_t starttime = clock();
		TR << "entering perturb steps" << std::endl;

		core::pose::Pose const saved_input_pose( pose ); //used to return sidechains later

		protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap(core::chemical::CENTROID);
		typeset_swap.apply( pose );
		//centroid score
		TR << "centroid score of starting PDB: " << (*centroid_scorefunction_)(pose) << std::endl;
		centroid_scorefunction_->show( TR, pose );
		TR << std::flush; //show doesn't flush the buffer

		/*
		perturbmover of either type
		minimize every once in a while
		MC evaluate
		*/

		////////////////////////////////////backbone_mover_cent/////////////////////////////////////
		protocols::moves::RandomMoverOP backbone_mover_cen( new protocols::moves::RandomMover() );

		using protocols::simple_moves::SmallMover;
		using protocols::simple_moves::BackboneMoverOP;
		protocols::simple_moves::BackboneMoverOP small_mover_cen = new protocols::simple_moves::SmallMover(movemap_, 0.8, 0);
		small_mover_cen->angle_max( 'H', 180.0 );
		small_mover_cen->angle_max( 'E', 180.0 );
		small_mover_cen->angle_max( 'L', 180.0 );
		backbone_mover_cen->add_mover(small_mover_cen, 1.0);

		protocols::simple_moves::BackboneMoverOP shear_mover_cen = new protocols::simple_moves::ShearMover(movemap_, 0.8, 0);
		shear_mover_cen->angle_max( 'H', 180.0 );
		shear_mover_cen->angle_max( 'E', 180.0 );
		shear_mover_cen->angle_max( 'L', 180.0 );
		//backbone_mover_cen->add_mover(shear_mover_cen, 1.0); //not yet

		if ( fragset3mer_ ) { //if we have fragments
			using protocols::simple_moves::ClassicFragmentMover;
			protocols::simple_moves::ClassicFragmentMoverOP frag_mover = new ClassicFragmentMover(fragset3mer_, movemap_);
			frag_mover->enable_end_bias_check(false);
			backbone_mover_cen->add_mover(frag_mover, 0.5);
		}

		utility::vector1< protocols::moves::MoverOP > loop_movers;
		if ( !loops_.empty() ) { //if we have loops to deal with
			for ( protocols::loops::Loops::const_iterator it(loops_.begin()), end(loops_.end()); it != loops_.end(); it++ ) {
				//make kinematic mover
				using protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP;
				using protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP;
				using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
				KinematicMoverOP kin_mover( new KinematicMover() );//temperature, default 0.8
				KinematicMoverCAP kin_mover_cap(*kin_mover);

				//set up kinematic mover - this is borrowed from loops_main.cc, perturb_one_loop_with_alc(), SVN 24219, #932-946
				kin_mover->set_vary_bondangles( false );
				kin_mover->set_sample_nonpivot_torsions( true );
				kin_mover->set_rama_check( true );
				kin_mover->set_idealize_loop_first( false );

				//make vicinity-sampling kinematic perturber
				using protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber;
				using protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP;
				VicinitySamplingKinematicPerturberOP TsamplingKP( new VicinitySamplingKinematicPerturber(kin_mover_cap));
				//TsamplingKP->set_sample_vicinity( true );
				TsamplingKP->set_degree_vicinity( 35 );
				kin_mover->set_perturber(TsamplingKP);

				using protocols::loops::loop_closure::kinematic_closure::KinematicWrapperOP;
				using protocols::loops::loop_closure::kinematic_closure::KinematicWrapper;
				KinematicWrapperOP kin_wrapper( new KinematicWrapper(kin_mover, *it));
				backbone_mover_cen->add_mover(kin_wrapper, 0.5);
				loop_movers.push_back(kin_wrapper);
				TR << "creating KinematicWrapper with loop " << *it << std::endl;

			}//for each loop
		}//if loops exist

		/////////////////////////minimizer mover/////////////////////////////////////////
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		protocols::simple_moves::MinMoverOP min_mover_cen = new protocols::simple_moves::MinMover(
			movemap_,
			centroid_scorefunction_,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true /*use_nblist*/ );

		/////////////////////////Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		using protocols::moves::MonteCarlo;
		using protocols::moves::MonteCarloOP;
		using basic::options::option;
		using namespace basic::options::OptionKeys::AnchoredDesign;
		MonteCarloOP mc_cen( new MonteCarlo( pose, *centroid_scorefunction_, option[ perturb_temp ].value() ) );

		///////////////////////////////////for loop///////////////////////////////////////////////////
		protocols::moves::PDBDumpMover cen_out("cen_cycle");
		core::Size const perturb_applies = option[ perturb_cycles ].value(); //default 5
		core::Size shear_on_cyc(core::Size(core::Real(perturb_applies) * option[ basic::options::OptionKeys::shear_on ]));
		if ( shear_on_cyc == 0 ) shear_on_cyc = 1; //0 should mean shear on immediately, but the if below never sees 0.
		TR << "shear on at " << shear_on_cyc << std::endl;
		TR << "   Current     Low    total cycles =" << perturb_applies << std::endl;
		for ( core::Size i = 1; i <= perturb_applies; ++i ) {
			if ( i == shear_on_cyc ) backbone_mover_cen->add_mover(shear_mover_cen, 1.0);
			if ( (i % 20 == 0) || (i == perturb_applies) ) min_mover_cen->apply(pose);
			else backbone_mover_cen->apply(pose);

			if ( option[ debug ] ) cen_out.apply(pose); //check trajectory
			mc_cen->boltzmann(pose);

			TR << i << "  " << mc_cen->last_accepted_score() << "  " << mc_cen->lowest_score() << std::endl;

			//constraint report
			//using core::scoring::atom_pair_constraint;
			//core::Real const cstscore_in(pose.energies().total_energies()[atom_pair_constraint] * pose.energies().weights()[atom_pair_constraint]);
			//TR << "cst score " << cstscore_in << std::endl;
		}//end the exciting for loop
		mc_cen->recover_low(pose);

		//filter based on constraints score - if not less than 1 (close to 0), cancel this trajectory
		// using core::scoring::atom_pair_constraint;
		// core::Real const cstscore(pose.energies().total_energies()[atom_pair_constraint] * pose.energies().weights()[atom_pair_constraint]);
		// if (cstscore > 1.0) {
		//  TR << "centroid constraints not satisfied; final constraint score: " << cstscore << ", restarting centroid" << std::endl;
		//  set_last_move_status(protocols::moves::FAIL_RETRY);
		//  return;
		// }

		//dump centroid-stage result pose
		if ( basic::options::option[basic::options::OptionKeys::AnchoredDesign::perturb_show ].value() ) {
			using namespace protocols::jd2;
			JobCOP job_me( JobDistributor::get_instance()->current_job() );
			JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "centroid");
		}

		//show centroid score (duplicates last line above)
		TR << "centroid score of final perturbed PDB: " << (*centroid_scorefunction_)(pose) << std::endl;
		centroid_scorefunction_->show( TR, pose );
		TR << std::flush; //show doesn't flush the buffer

		clock_t stoptime = clock();
		TR << "One perturb took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;
		TR << "perturb steps complete" << std::endl;
		starttime = clock();

		///////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////fullatom///////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////

		protocols::simple_moves::ReturnSidechainMover return_sidechains( saved_input_pose );
		return_sidechains.apply( pose );

		//remove centroid constraints; add fullatom constraints
		pose.remove_constraints();
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose); //protected internally if no csts

		/////////////////////////////generate full repack&minimize mover//////////////////////////////
		protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
		pack_mover->task_factory( task_factory_ );
		pack_mover->score_function( fullatom_scorefunction_ );

		protocols::simple_moves::MinMoverOP min_mover_fa = new protocols::simple_moves::MinMover(
			movemap_,
			fullatom_scorefunction_,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true /*use_nblist*/ );

		//definitely want sidechain minimization here
		using protocols::simple_moves::TaskAwareMinMoverOP;
		using protocols::simple_moves::TaskAwareMinMover;
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover_fa = new protocols::simple_moves::TaskAwareMinMover(min_mover_fa, task_factory_);

		/////////////////////////repack/minimize once to fix sidechains//////////////////////////////////
		// TR << "packing" << std::endl;
		pack_mover->apply(pose);
		// TR << "minimizing" << std::endl;
		TAmin_mover_fa->apply(pose);

		//////////////////////////////////////// backbone mover/////////////////////////////////////////
		protocols::moves::RandomMoverOP backbone_mover_fa( new protocols::moves::RandomMover() );

		protocols::simple_moves::BackboneMoverOP small_mover_fa = new protocols::simple_moves::SmallMover(movemap_lesstail_, 0.8, 0);
		small_mover_fa->angle_max( 'H', 4.0 );
		small_mover_fa->angle_max( 'E', 4.0 );
		small_mover_fa->angle_max( 'L', 4.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover_fa = new protocols::simple_moves::ShearMover(movemap_lesstail_, 0.8, 0);
		shear_mover_fa->angle_max( 'H', 4.0 );
		shear_mover_fa->angle_max( 'E', 4.0 );
		shear_mover_fa->angle_max( 'L', 4.0 );

		backbone_mover_fa->add_mover(small_mover_fa, 1.0);
		backbone_mover_fa->add_mover(shear_mover_fa, 1.0);

		/////////////////////loops as needed///////////////////////////////
		if ( !loops_.empty() && !loop_movers.empty() ) { //if we have loops to deal with
			for ( int i(1); i<=loop_movers.size(); ++i ) {
				backbone_mover_fa->add_mover(loop_movers[i], 0.5);  //same loop movers from centroid mode should work fine
			}
		}

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		MonteCarloOP mc_fa( new MonteCarlo( pose, *fullatom_scorefunction_, option[ refine_temp ].value() ) );

		/////////////////////////////////rotamer trials mover///////////////////////////////////////////
		using protocols::simple_moves::RotamerTrialsMoverOP;
		using protocols::simple_moves::EnergyCutRotamerTrialsMover;
		protocols::simple_moves::RotamerTrialsMoverOP rt_mover(new protocols::simple_moves::EnergyCutRotamerTrialsMover(
			fullatom_scorefunction_,
			task_factory_,
			mc_fa,
			0.01 /*energycut*/ ) );

		/////////////////////////////////////////refine loop///////////////////////////////////////////
		core::Size const refine_applies = option[ refine_cycles ].value(); //default 5
		core::Size const repack_cycles = option[ refine_repack_cycles ].value();
		core::Size const min_cycles = repack_cycles/2;
		core::Size const switch_movemaps(core::Size(core::Real(refine_applies) * option[ basic::options::OptionKeys::short_tail_off ]));
		TR << "   Current     Low    total cycles =" << refine_applies << std::endl;
		for ( core::Size i(1); i <= refine_applies; ++i ) {
			if ( i == switch_movemaps ) {
				small_mover_fa->movemap(movemap_);
				shear_mover_fa->movemap(movemap_);
			}
			if ( (i % repack_cycles == 0) || (i == refine_applies) ) { //full repack
				pack_mover->apply(pose);
				TAmin_mover_fa->apply(pose);
			} else if ( i % min_cycles == 0 ) { //minimize
				TAmin_mover_fa->apply(pose);
			} else {
				backbone_mover_fa->apply(pose);
				rt_mover->apply(pose);
			}

			mc_fa->boltzmann(pose);
			TR << i << "  " << mc_fa->last_accepted_score() << "  " << mc_fa->lowest_score() << std::endl;
		}//end the exciting for loop
		mc_fa->recover_low( pose );

		//let's store some energies/etc of interest
		//this code is specific to the E2/RING/E3 system for which this code was written; it is refactored elsewhere
		// BARAK: this line should be commented out if applied to other systems
		//SML 2/1/11: it's now under commandline control
		if ( basic::options::option[ basic::options::OptionKeys::publication].value() ) create_extra_output(pose, fullatom_scorefunction_);

		(*fullatom_scorefunction_)(pose);
		set_last_move_status(protocols::moves::MS_SUCCESS); //this call is unnecessary but let's be safe
		return;
	}

	virtual std::string get_name() const { return "FloppyTailMover"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new FloppyTailMover;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return true; }

private:
	core::Size start_;
	core::Size stop_;
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP centroid_scorefunction_;
	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP movemap_;
	core::kinematics::MoveMapOP movemap_lesstail_;
	/// @brief stored so that it can be generated in the init_on_new_input function
	core::kinematics::FoldTreeOP foldtree_;
	core::fragment::ConstantLengthFragSetOP fragset3mer_;

	protocols::loops::Loops loops_;
};

typedef utility::pointer::owning_ptr< FloppyTailMover > FloppyTailMoverOP;

int main( int argc, char* argv[] )
{

	try {

		using basic::options::option;
		using namespace basic::options::OptionKeys;
		option.add( flexible_start_resnum, "start res for flexible region").def(180);
		option.add( flexible_stop_resnum, "stop res for flexible region").def(0);
		option.add( flexible_chain, "chain for flexible region").def("C");
		option.add( shear_on, "fraction of perturb moves when shear turns on (0.5 = halfway through").def(1.0/3.0);
		option.add( short_tail_fraction, "what fraction of the flexible segment is used in the short-tail section of refinement (not compatible with non-terminal flexible regions)").def(1.0/3.0);
		option.add( short_tail_off, "fraction of refine cycles where movemap reverts to full tail (0.5 = halfway through)").def(0.1);
		option.add( pair_off, "turn off Epair electrostatics term").def(false);
		option.add( publication, "output statistics used in publication.  TURN OFF if not running publication demo.").def(false);
		option.add( C_root, "Reroot the fold_tree to the C-terminus.  If your flexible region is N-terminal, or closer to the first half of the pose, this will speed computation.").def(false);

		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(new FloppyTailMover);

		TR << "************************d**o**n**e**************************************" << std::endl;

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
