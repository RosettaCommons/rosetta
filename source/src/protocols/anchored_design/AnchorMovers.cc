// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/anchored_design/AnchorMovers.cc
/// @brief Anchored design movers (protocol-level, perturb, refine)
/// @author Steven Lewis smlewi@gmail.com

// Unit Headers
#include <protocols/anchored_design/AnchorMovers.hh>

// Package Headers
#include <protocols/anchored_design/AnchorMoversData.hh>
#include <protocols/analysis/LoopAnalyzerMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/loops/Loops.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/fragment/FragSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

//needed for a benchmarking thing
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <basic/MetricValue.hh>

//movers and accessories
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/moves/Mover.fwd.hh> //MoverOP
#include <protocols/simple_moves/BackboneMover.hh> //SmallMover
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> //typeset swapping
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <ObjexxFCL/FArray1D.hh> //necessary for fold tree tricks, RMS - AUGH!
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/exit.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <iostream>
#include <sstream>
#include <list>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.io.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer T_design( "protocols.AnchoredDesign.AnchoredDesignMover" );
static THREAD_LOCAL basic::Tracer T_perturb( "protocols.AnchoredDesign.AnchoredPerturbMover" );
static THREAD_LOCAL basic::Tracer T_refine( "protocols.AnchoredDesign.AnchoredRefineMover" );
static THREAD_LOCAL basic::Tracer T_shared( "protocols.AnchoredDesign.Anchor_Movers" );

namespace protocols {
namespace anchored_design {

std::string const EMPTY_STRING("");
int const ANCHOR_TARGET(1); //jump between anchor and target

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////local helper functions - generally these only run when the debug flag is passed/////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief a local helper function, not in the header
void debug_dump_pose(
	core::pose::Pose const & pose,
	core::pose::PoseOP posecopy,
	std::string const & tag,
	protocols::moves::MoverOP mover)
{
	using namespace protocols::jd2;
	mover->apply(*posecopy);
	JobDistributor::get_instance()->job_outputter()->other_pose(JobDistributor::get_instance()->current_job(), *posecopy, tag);
	*posecopy = pose;
}

/// @brief a local helper function for debugging, not in the header
void dump_cutpoint_info( core::pose::Pose const & pose)
{
	core::Size nres(pose.total_residue());
	T_shared << "pose says 'is_cutpoint(#)': ";
	for ( core::Size i(1); i <= nres; ++i ) {
		if ( pose.fold_tree().is_cutpoint( i ) ) T_shared << i << " ";
	}
	T_shared << std::endl;

	T_shared << "pose says 'has_variant_type CUTPOINT_LOWER': ";
	for ( core::Size i(1); i <= nres; ++i ) {
		if ( pose.residue(i).has_variant_type( core::chemical::CUTPOINT_LOWER) ) T_shared << i << " ";
	}
	T_shared << std::endl;

	T_shared << "pose says 'has_variant_type CUTPOINT_UPPER': ";
	for ( core::Size i(1); i <= nres; ++i ) {
		if ( pose.residue(i).has_variant_type( core::chemical::CUTPOINT_UPPER) ) T_shared << i << " ";
	}
	T_shared << std::endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////AnchoredDesignMover////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//@brief constructor with arguments
AnchoredDesignMover::AnchoredDesignMover( protocols::anchored_design::AnchorMoversDataOP interface_in ) :
	Mover(),
	interface_(std::move( interface_in )),
	RMSD_only_this_pose_( /* 0 */ ),
	IAM_( protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(ANCHOR_TARGET) ) ),
	rmsd_(false),
	RMSD_only_this_(EMPTY_STRING),
	delete_interface_native_sidechains_(false),
	show_extended_(false),
	randomize_input_sequence_(false),
	vary_cutpoints_(false),
	refine_only_(false),
	filter_score_(std::numeric_limits<core::Real>::max()),
	use_filter_score_(false),
	filter_SASA_(0.0),
	use_filter_SASA_(false),
	use_filter_omega_(false),
	autoinitialize_(true),
	init_for_input_yet_(false)
{
	protocols::moves::Mover::type( "AnchoredDesign" );
}

//@brief constructor with no arguments
AnchoredDesignMover::AnchoredDesignMover() :
	interface_( /* 0 */ ), //NULL pointer
	RMSD_only_this_pose_( /* 0 */ ),
	IAM_( protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover(ANCHOR_TARGET) ) ),
	rmsd_(false),
	RMSD_only_this_(EMPTY_STRING),
	delete_interface_native_sidechains_(false),
	show_extended_(false),
	randomize_input_sequence_(false),
	vary_cutpoints_(false),
	refine_only_(false),
	filter_score_(std::numeric_limits<core::Real>::max()),
	use_filter_score_(false),
	filter_SASA_(0.0),
	use_filter_SASA_(false),
	use_filter_omega_(false),
	autoinitialize_(true),
	init_for_input_yet_(false)
{
	protocols::moves::Mover::type( "AnchoredDesign" );
}

/// @brief copy ctor
AnchoredDesignMover::AnchoredDesignMover( AnchoredDesignMover const & rhs ) :
	Mover(rhs)
{
	*this = rhs;
}

/// @brief assignment operator
AnchoredDesignMover & AnchoredDesignMover::operator=( AnchoredDesignMover const & rhs ){

	//abort self-assignment
	if ( this == &rhs ) return *this;

	interface_ = rhs.interface_->clone();
	//this is a shallow copy because it is only supposed to have one value anyway (only this!), what does it mean if you start mutating it?  I don't know
	RMSD_only_this_pose_ = rhs.RMSD_only_this_pose_;
	//this is a shallow copy because the Mover can't be modified anyway, and it doesn't have a copy ctor at the moment, and it's not going to get one in the near future
	IAM_ = rhs.IAM_;
	rmsd_ = rhs.get_rmsd();
	RMSD_only_this_ = rhs.get_RMSD_only_this();
	delete_interface_native_sidechains_ = rhs.get_delete_interface_native_sidechains();
	show_extended_ = rhs.get_show_extended();
	randomize_input_sequence_ = rhs.get_randomize_input_sequence();
	vary_cutpoints_ = rhs.get_vary_cutpoints();
	refine_only_ = rhs.get_refine_only();
	filter_score_ = rhs.get_filter_score();
	use_filter_score_ = rhs.use_filter_score_;
	filter_SASA_ = rhs.get_filter_SASA();
	use_filter_SASA_ = rhs.use_filter_SASA_;
	use_filter_omega_ = rhs.get_filter_omega();
	autoinitialize_ = rhs.get_autoinitialize();
	init_for_input_yet_ = rhs.init_for_input_yet_;
	return *this;
}

AnchoredDesignMover::~AnchoredDesignMover() = default;

void AnchoredDesignMover::init_on_new_input(core::pose::Pose const & pose) {
	//don't run this function twice
	init_for_input_yet_ = true;

	//If the interface_ object doesn't exist yet, we must create one, giving it a pose to help it initialize
	//it will exist if this object was created via the AnchorMoversData-supplying constructor
	if ( !interface_ ) {
		interface_ = protocols::anchored_design::AnchorMoversDataOP( new protocols::anchored_design::AnchorMoversData(pose) );
	}

	//If nobody told us not to autoinitialize, read the options system for our other data
	if ( autoinitialize_ ) read_options();

	//If we are in RMSD_only_this_ mode, we will need to set up that cached comparison pose
	if ( RMSD_only_this_ != EMPTY_STRING ) {
		core::pose::Pose dummy;
		core::import_pose::pose_from_file(dummy, RMSD_only_this_, core::import_pose::PDB_file);
		RMSD_only_this_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(dummy) ) );
	}

	return;
}

/// @details AnchoredDesignMover is mostly a container for other movers for the anchored design protocol.
void AnchoredDesignMover::apply( core::pose::Pose & pose )
{
	clock_t starttime = clock();

	if ( !init_for_input_yet_ ) init_on_new_input(pose);

	core::pose::PoseCOP start_pose(nullptr);

	//pre-pre-processing
	if ( rmsd_ ) {
		start_pose = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose) ) );
	}

	if ( RMSD_only_this_ != EMPTY_STRING ) {
		core::pose::Pose dummy;
		core::import_pose::pose_from_file(dummy, RMSD_only_this_, core::import_pose::PDB_file);
		start_pose = RMSD_only_this_pose_;
	} else { //if RMSD_only_this is active, we skip all the meat steps


		//pre-processing
		interface_->pick_new_cutpoints(vary_cutpoints_); //verifies that cutpoints are legal w/r/t anchor; chooses new ones if not; picks new random ones if bool is true
		set_fold_tree_and_cutpoints( pose );
		forget_initial_loops( pose ); //checks loop object is_extended internally

		if ( randomize_input_sequence_ ) randomize_input_sequence(pose);
		if ( delete_interface_native_sidechains_ ) delete_interface_native_sidechains(pose);

		if ( interface_->get_anchor_noise_constraints_mode() ) {
			interface_->anchor_noise_constraints_setup(pose);
			perturb_anchor(pose);
		}

		if ( show_extended_
				|| randomize_input_sequence_
				|| delete_interface_native_sidechains_
				|| interface_->get_anchor_noise_constraints_mode() ) {
			using namespace protocols::jd2;
			JobDistributor::get_instance()->job_outputter()->other_pose(JobDistributor::get_instance()->current_job(), pose, "preprocessed");
		}

		//processing
		if ( !refine_only_ ) {
			protocols::anchored_design::AnchoredPerturbMover anchor_perturb( interface_ );
			anchor_perturb.apply( pose );
		}

		//do not want constraints in centroid mode - for one, they're probably fullatom, and for two, the centroid scorefunction isn't so hot anyway
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose ); //protected internally if no constraints

		protocols::anchored_design::AnchoredRefineMover anchor_refine( interface_ );
		anchor_refine.apply( pose );

		//post-processing
		protocols::analysis::LoopAnalyzerMover LAM( interface_->loops());
		LAM.apply( pose );

		//load sequence into Job output
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string("SEQUENCE: " + pose.sequence());

		//load interface analysis info
		IAM_->apply( pose );

		//load unconstrainted score
		if ( pose.constraint_set()->has_constraints() ) {
			pose.constraint_set()->show_definition(T_design, pose);
			core::pose::Pose nocstcopy(pose);
			nocstcopy.remove_constraints();
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("no_cst_total", ((*interface_->get_fullatom_scorefunction())(nocstcopy)));
		}

		clock_t stoptime = clock();
		T_design << "One perturb/refine took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;
	} //end RMSD_only_this block

	//post-post-processing
	calculate_rmsd(pose, start_pose);

	filter(pose);
	(*interface_->get_fullatom_scorefunction())(pose);
	return;
}//AnchoredDesignMover::apply

void
AnchoredDesignMover::calculate_rmsd( core::pose::Pose const & pose, core::pose::PoseCOP start_pose ){

	if ( rmsd_ ) {
		core::Real const rmsd(core::scoring::CA_rmsd(pose, *start_pose));
		T_design << "CA_sup_RMSD for this trajectory: " << rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("CA_sup_RMSD", rmsd);

		std::list<core::Size> loops_for_rmsd;
		for ( core::Size i(1), num_loops(interface_->num_loops()); i <= num_loops; ++i ) {
			core::Size const loopstart(interface_->loop(i).start()), loopend(interface_->loop(i).stop());
			for ( core::Size j = loopstart; j <= loopend; ++j ) {
				loops_for_rmsd.insert(loops_for_rmsd.end(), j);
			}
		}

		//   T_design << "loop residues for RMSD:";
		//   for( std::list<core::Size>::const_iterator it(loops_for_rmsd.begin()), end(loops_for_rmsd.end()); it != end; ++it){
		//    T_design << " " << *it;
		//   }
		//   T_design << std::endl;

		core::Real const loop_rmsd(core::scoring::CA_rmsd(pose, *start_pose, loops_for_rmsd));
		T_design << "loop_CA_sup_RMSD for this trajectory: " << loop_rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("loop_CA_sup_RMSD", loop_rmsd);

		//chain 1 RMSD
		core::Real const ch1_sup_rmsd(core::scoring::CA_rmsd(pose, *start_pose, pose.conformation().chain_begin(1), pose.conformation().chain_end(1)));
		T_design << "chain 1 sup_RMSD for this trajectory: " << ch1_sup_rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("ch1_CA_sup_RMSD", ch1_sup_rmsd);

		//chain 2 RMSD
		core::Real const ch2_sup_rmsd(core::scoring::CA_rmsd(pose, *start_pose, pose.conformation().chain_begin(2), pose.conformation().chain_end(2)));
		T_design << "chain 2 sup_RMSD for this trajectory: " << ch2_sup_rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("ch2_CA_sup_RMSD", ch2_sup_rmsd);

		//also need no-superimpose RMSDs for chain1/2 to get lever arm effect
		ObjexxFCL::FArray1D_bool ch1(pose.total_residue(), false), ch2(pose.total_residue(), false);

		for ( core::Size i(pose.conformation().chain_begin(1)), end(pose.conformation().chain_end(1)); i<=end; ++i ) ch1(i)=true;
		//chain 1 RMSD
		core::Real const ch1_rmsd(core::scoring::rmsd_no_super_subset(pose, *start_pose, ch1, core::scoring::is_protein_CA));
		T_design << "chain 1 RMSD for this trajectory: " << ch1_rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("ch1_CA_RMSD", ch1_rmsd);

		for ( core::Size i(pose.conformation().chain_begin(2)), end(pose.conformation().chain_end(2)); i<=end; ++i ) ch2(i)=true;
		//chain 2 RMSD
		core::Real const ch2_rmsd(core::scoring::rmsd_no_super_subset(pose, *start_pose, ch2, core::scoring::is_protein_CA));
		T_design << "chain 2 RMSD for this trajectory: " << ch2_rmsd << std::endl;
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("ch2_CA_RMSD", ch2_rmsd);

		//This next section is meant purely for benchmarking purposes.  It assumes the moving chain is chain 2, and that the interface-definition calculator created by the standard AnchorClass default TaskFactory exists.  It will use rmsd_with_super_subset on interface backbone atoms, because that's what docking does.
		//chain 2 interface RMSD
		if ( true ) {

			//get interface set from calculator
			std::string const & interface_calc(interface_->interface_calc());
			//find the set of residues
			typedef std::set< core::Size > SizeSet;
			basic::MetricValue< SizeSet > mv_sizeset;
			start_pose->metric(interface_calc, "interface_residues", mv_sizeset);
			SizeSet const & sizeset(mv_sizeset.value());

			//convert this set into the type needed by rmsd_with_super_subset
			T_design << "interface for rmsd";
			ObjexxFCL::FArray1D_bool is_interface( pose.total_residue(), false );
			for ( unsigned long it : sizeset ) {
				is_interface(it) = true;
				T_design << " " << it;
			}
			T_design << std::endl;

			//ready to calculate
			core::Real const I_sup_bb_rmsd(core::scoring::rmsd_with_super_subset(pose, *start_pose, is_interface, core::scoring::is_protein_backbone));
			T_design << "interface backbone RMSD for this trajectory: " << I_sup_bb_rmsd << std::endl;
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair("I_sup_bb_RMSD", I_sup_bb_rmsd);

		}

	}
	return;
}

/// @details This is crazy.  Sometimes loop redesign is overly biased by the starting loop sequence, because the centroid phase won't change the sequence, and the rama term may remember the conformation to well due to sequence preferences.  This restricts sampling.  So, this function creates a random sequence in the designable positions, based on what the design PackerTask thinks they can become.  This code special cases histidine (to avoid double-allowing histidine, based on its two protonations).
void AnchoredDesignMover::randomize_input_sequence(core::pose::Pose & pose ) const {

	T_design << "entering randomize_input_sequence" << std::endl;

	//Get copy of "right" PackerTask
	core::pack::task::PackerTaskCOP oldtask(interface_->get_task(pose));

	//there isn't much we can do to use the old task due to accumulation of status; also (and more importantly) using a task to pack is ruined by uneven numbers of rotamers (arginine would be favored over glycine).

	core::Size const nres(pose.total_residue());
	for ( core::Size i(1); i<=nres; ++i ) { //for all positions
		//if this position is designable, do something
		using core::pack::task::ResidueLevelTask;
		ResidueLevelTask const & rlt(oldtask->residue_task(i));
		if ( rlt.being_designed() ) {
			//get the list
			ResidueLevelTask::ResidueTypeCOPList types(rlt.allowed_residue_types());

			//print possibilities before histidine check
			T_design << "before HIS/D check, position " << i;

			for ( auto & type : types ) {
				T_design << " " << type->name();// << std::endl;
			}
			T_design << std::endl;

			//sweep for histidines
			core::Size num_histidines(0); //count histidines we find; if more than 2 explode
			core::chemical::ResidueTypeCOP histidine; //store the most recent histidine we find as we go
			for ( auto
					allowed_iter = types.begin(),
					iter_next = types.begin(),
					allowed_end = types.end();
					allowed_iter != allowed_end;  /* no increment: deletion + iterator incrementing = segfault! */ ) {
				iter_next = allowed_iter;
				++iter_next;

				if ( (*allowed_iter)->aa() == core::chemical::aa_his ) { //we only want to look at histidines
					++num_histidines;
					histidine=*allowed_iter;
					types.erase( allowed_iter );
				}
				allowed_iter = iter_next;
			}//histidine removal scan
			if ( num_histidines > 0 && num_histidines < 3 ) {
				types.push_back(histidine); //if we removed 1 or 2
			} else if ( num_histidines < 3 ) {
				utility_exit_with_message("removed more than 2 histidines in AnchoredDesign::randomize_input_sequence, something is wrong");
			}

			//print possibilities after histidine check
			T_design << "after HIS/D check, position " << i;

			for ( auto & type : types ) {
				T_design << " " << type->name();
			}
			T_design << std::endl;

			//now that that's out of the way, pick a ResidueType at random
			core::Size const num_types(types.size());
			core::Size const chosen_type_index(numeric::random::rg().random_range(1, num_types));
			auto iter = types.begin();
			for ( core::Size add(1); add<chosen_type_index; ++add ) ++iter;
			core::chemical::ResidueTypeCOP chosen_type(*iter);

			T_design << "chose at position " << i << chosen_type->name() << std::endl;
			T_design << "chose at position " << i << chosen_type->name() << std::endl;

			//do a replace_replace with the chosen type
			using namespace core::conformation;
			ResidueCOP new_residue(ResidueFactory::create_residue( *chosen_type, pose.residue(i), pose.conformation()));
			pose.replace_residue(i, *new_residue, true);

		}//if being designed
		//if not being designed, do nothing
	}//for all residues

	return;
}

/// @details For benchmarking, it is a minor sin to allow native sidechains to leak through from the starting structure.  AnchoredDesign runs best with use_input_sc because it does sidechain minimization, not because it needs starting sidechains.  This function deletes the native sidechains by repacking the interface with use_input_sc forcibly off.  DO NOT USE THIS FUNCTION for proper designs - it is meant for a specific benchmarking purpose.
void AnchoredDesignMover::delete_interface_native_sidechains(core::pose::Pose & pose ) const {

	//create a TaskFactory for the deletion.  Notice that it DOES NOT read all user inputs (ignoring the resfile and command line) and thus can only be used when you satisfy its assumptions.  The assumptions are that you are A) doing fixed-sequence benchmarking, and B) starting from the correct interface structure.  The factory is set up to detect the interface and repack those side chains, excepting the anchor, to with many rotamers.  It does NOT pay attention to the flags that allow repacking of the anchor.  It ASSUMES that the AnchorMoversData class has pregenerated some PoseMetricCalculators.
	using namespace core::pack::task;
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );

	//some code copied from AnchorClass.cc:AnchorMoversData::set_unset_packertask_factory.  I chose to copy (gasp) because I do want the two to be uncoupled; if the other code evolves this should stay the same.  Of course, the PoseMetricCalculators are still a dependency...
	//operation to detect interface/loops; depends on preexistence of these calculators
	std::string const & interface_calc(interface_->interface_calc());
	std::string const & neighborhood_calc(interface_->neighborhood_calc());
	utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
	calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
	calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

	using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
	operation::TaskOperationOP rbcop( new RestrictByCalculatorsOperation( calcs_and_calcns ) );

	tf->push_back(rbcop);

	//operation to protect anchor
	operation::PreventRepackingOP prop( new operation::PreventRepacking );
	for ( core::Size i(interface_->anchor_start()); i<= interface_->anchor_end(); ++i ) {
		prop->include_residue(i);
	}

	tf->push_back(prop);

	//operations to allow ex flags to saturation
	for ( core::Size i(1), end(pose.total_residue()); i<=end; ++i ) {
		tf->push_back(operation::RotamerExplosionOP( new core::pack::task::operation::RotamerExplosion(i, core::pack::task::EX_ONE_STDDEV, 4) ));
		tf->push_back(operation::ExtraChiCutoffOP( new core::pack::task::operation::ExtraChiCutoff(i, 0) ));
		//tf->push_back(operation::IncludeCurrentOP( new core::pack::task::operation::IncludeCurrent())); //this is exactly what we do not want - useful for testing that it worked right...
	}

	//we DO NOT WANT design
	tf->push_back(operation::TaskOperationCOP( new operation::RestrictToRepacking() ));

	//print a copy of the task for double checking
	//T_shared << *(tf->create_task_and_apply_taskoperations(pose)) << std::endl;

	//pose.dump_pdb("pre_delete_interface_sidechains_test.pdb");
	//create a PackRotamersMover and do it
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
	pack_mover->task_factory( tf );
	pack_mover->score_function( interface_->get_fullatom_scorefunction() );
	pack_mover->apply(pose);
	//pose.dump_pdb("post_delete_interface_sidechains_test.pdb");
	return;
}

/// @details when using anchor_noise_constraints_mode, this function perturbs the initial anchor a bit.  Requested test by a reviewer.
void AnchoredDesignMover::perturb_anchor( core::pose::Pose & pose ) const {

	//I think if we take the Stubs out of the existing Jump, all we have to do is change the vector part
	using core::kinematics::Stub;
	using core::kinematics::Jump;
	Stub const & S1( pose.conformation().upstream_jump_stub( ANCHOR_TARGET ) );
	Stub const & S2( pose.conformation().downstream_jump_stub( ANCHOR_TARGET ) );
	core::Vector const & translation(S2.v);

	//store xyz to check momentarily
	core::Size const anchor(interface_->anchor_start());
	core::Size const CA(pose.residue_type(anchor).atom_index("CA"));
	core::id::AtomID const anchor_ID(core::id::AtomID(CA, anchor));
	core::Vector const original_anchor_xyz(pose.xyz(anchor_ID));

	//We want to make a change of plus or minus one angstrom.  numeric::random::rg().uniform returns a range of 0 to 1, so the formula 2N-1 converts it to a range of -1 to 1 perturbation - we add this to the original coordinate.
	core::Real const new_x(translation.x() + ((numeric::random::rg().uniform()*2.0) - 1.0));
	core::Real const new_y(translation.y() + ((numeric::random::rg().uniform()*2.0) - 1.0));
	core::Real const new_z(translation.z() + ((numeric::random::rg().uniform()*2.0) - 1.0));

	T_design << "perturb_anchor: old/new, x->y->z\n"
		<< translation.x() << " " << new_x << "\n"
		<< translation.y() << " " << new_y << "\n"
		<< translation.z() << " " << new_z << std::endl;

	core::Vector const new_trans(new_x, new_y, new_z);

	pose.set_jump(ANCHOR_TARGET, Jump(S1, Stub( S2.M, new_trans ) ));

	core::Vector const new_anchor_xyz(pose.xyz(anchor_ID));

	T_design << "perturb_anchor: old and new CA positions:\n"
		<< original_anchor_xyz << '\n'
		<< new_anchor_xyz << std::endl;

	return;
}

std::string
AnchoredDesignMover::get_name() const {
	return "AnchoredDesignMover";
}

protocols::moves::MoverOP AnchoredDesignMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AnchoredDesignMover );
}

protocols::moves::MoverOP AnchoredDesignMover::clone() const {
	return protocols::moves::MoverOP( new AnchoredDesignMover(*this) );
}

bool AnchoredDesignMover::reinitialize_for_each_job() const { return false; }

/// @details generally returns true; will return false for RMSD_only_this mode
bool AnchoredDesignMover::reinitialize_for_new_input() const { return (RMSD_only_this_ == EMPTY_STRING); }

void AnchoredDesignMover::filter( core::pose::Pose & pose ){

	std::ostringstream failure;
	bool fail(false);

	if ( use_filter_score_ ) {
		core::Real const pscore((*interface_->get_fullatom_scorefunction())(pose));
		if ( pscore > filter_score_ ) {
			failure << "failed total score filter; score " << pscore;
			fail = true;
		}
	}

	if ( use_filter_SASA_ && !fail ) {
		basic::MetricValue< core::Real > mv_delta_sasa;
		pose.metric("InterfaceSasaDefinition_1", "delta_sasa", mv_delta_sasa); //magic string: this calculator was created by the InterfaceAnalyzerMover
		if ( mv_delta_sasa.value() < filter_SASA_ ) {
			failure << "failed interface SASA filter; sasa " << mv_delta_sasa.value();
			fail = true;
		}
	}

	if ( use_filter_omega_ && !fail ) {
		core::Size const num_loops = interface_->num_loops();
		for ( core::Size i(1); i <= num_loops; ++i ) {
			core::Size const loopstart(interface_->loop(i).start()), loopend(interface_->loop(i).stop());
			for ( core::Size j = loopstart; j <= loopend; ++j ) {
				core::Real const omega(numeric::principal_angle_degrees(pose.omega(j)));
				if ( pose.residue(j).is_upper_terminus() ) continue; //omega will fail for c-term
				if ( (omega < 160) && (omega > -160) ) {
					failure << "failed omega at " << j << " " << omega;
					fail = true;
					i = num_loops + 1; //break out of parent loop
					break; //and this loop
				} //if omega is bad
			} //for residues in loop
		} //for each loop
	}//if filter active

	if ( !fail ) { //succeed case
		set_last_move_status(protocols::moves::MS_SUCCESS); //this call is unnecessary but let's be safe
		return;
	}

	//print failure
	T_design << failure.str() << std::endl;
	set_last_move_status(protocols::moves::FAIL_RETRY);
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details uses fold_tree::tree_from_jumps_and_cuts to make the necessary fold tree and apply it to the pose
void AnchoredDesignMover::set_fold_tree_and_cutpoints( core::pose::Pose & pose )
{
	core::Size const nres = pose.total_residue();

	//we'll need to know where the chainbreak is (we are assuming there is only one)
	core::Size const chain1end(pose.conformation().chain_end(1));
	if ( pose.conformation().num_chains() != 2 ) { //2 is the value we expect
		Warning() << "AnchoredDesign only tested with two chains; later chains may(?) be rigidly attached to chain 2 COOH, or not move at all, or who knows what";
	}

	core::Size const num_loops(interface_->num_loops());
	core::Size num_jumps(num_loops+1);//one extra for the inter-rigid-body jump

	//these may need to be redimensioned
	utility::vector1<int> jump_starts(num_jumps, 0), jump_stops(num_jumps, 0), cuts(num_jumps, 0);

	//using core::Real;
	core::Real anchorstart(interface_->anchor_start()), anchorend(interface_->anchor_end());
	core::Size anchormid(core::Size(std::ceil(anchorstart + (anchorend - anchorstart)/2.0)));

	Size anchor_jump_end( anchormid > chain1end ? chain1end : chain1end+1);
	jump_starts[ANCHOR_TARGET] = anchormid < anchor_jump_end ? anchormid : anchor_jump_end; //max
	jump_stops[ ANCHOR_TARGET] = anchormid > anchor_jump_end ? anchormid : anchor_jump_end; //min
	cuts[ANCHOR_TARGET]        = chain1end;
	assert(anchormid != anchor_jump_end);

	//  jumps(1, ANCHOR_TARGET) = anchormid < anchor_jump_end ? anchormid : anchor_jump_end; //max
	//  jumps(2, ANCHOR_TARGET) = anchormid > anchor_jump_end ? anchormid : anchor_jump_end; //min
	//  assert(anchormid != anchor_jump_end);
	//  cuts(ANCHOR_TARGET) = chain1end;

	int jump_num = ANCHOR_TARGET; //this reserves earlier jump numbers for the anchor
	//iterate over all loops
	for ( core::Size i(1); i <= num_loops; ++i ) {
		//if this is a terminal loop, don't mess with fold_tree, we won't need space in the vectors
		if ( pose.residue(interface_->loop(i).start()).is_terminus() //or
				|| pose.residue(interface_->loop(i).stop() ).is_terminus() ) {
			jump_starts.pop_back();
			jump_stops.pop_back();
			cuts.pop_back();
			--num_jumps; //one fewer jump to set up
			continue;
		}

		jump_num++;
		jump_starts[jump_num] = interface_->loop(i).start() - 1;
		jump_stops[ jump_num] = interface_->loop(i).stop()  + 1;
		cuts[jump_num] = interface_->loop(i).cut();

		//label cutpoints properly
		using core::pose::add_variant_type_to_pose_residue;
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, interface_->loop(i).cut() );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, interface_->loop(i).cut()+1 );
	}//over all loops

	//convert the vectors into FArrays
	assert(cuts.size() == num_jumps);
	ObjexxFCL::FArray2D<int> Fjumps(2, num_jumps);
	ObjexxFCL::FArray1D<int> Fcuts(num_jumps);
	for ( core::Size i(1); i<=num_jumps; ++i ) {
		Fjumps(1, i) = jump_starts[i];
		Fjumps(2, i) = jump_stops[i];
		Fcuts(i) = cuts[i];
	}

	//instantiate, create, output, apply our fold tree
	core::kinematics::FoldTree foldtree(nres);
	foldtree.tree_from_jumps_and_cuts(nres, num_jumps, Fjumps, Fcuts, 1, true);
	//foldtree.reorder(1);
	T_design << "anchored_design_fold_tree: " << foldtree << std::endl;

	//do it harder, better, faster, stronger
	core::kinematics::Edge const & anchoredge(foldtree.jump_edge(ANCHOR_TARGET));
	core::kinematics::Edge const new_anchoredge( anchoredge.start(), anchoredge.stop(), ANCHOR_TARGET, "CA", "CA", true);
	foldtree.delete_edge(anchoredge);
	foldtree.add_edge(new_anchoredge);
	foldtree.reorder(1);
	T_design << "anchored_design_fold_tree: " << foldtree << std::endl;

	pose.fold_tree(foldtree);
	dump_cutpoint_info(pose);
}//AnchoredDesignMover::set_fold_tree_and_cutpoints

/// @details implements the is_extended boolean for the Loop class.  If the boolean is on, this function resets phi/psi for those regions to be extended (-150phi, 150psi).  This function does not affect omega, or residues outside of defined loops, or the anchor.
void AnchoredDesignMover::forget_initial_loops( core::pose::Pose & pose ){

	//for each loop
	core::Size const num_loops(interface_->num_loops());
	for ( core::Size loop(1); loop <= num_loops; ++loop ) {
		//if loop has extended boolean set
		if ( interface_->loop(loop).is_extended() ) {
			//iterate through loop residues
			core::Size const loop_start(interface_->loop(loop).start()), loop_end(interface_->loop(loop).stop());
			for ( core::Size res(loop_start); res<=loop_end; ++res ) {
				//check movemap before alteration
				using namespace core::id;
				if ( interface_->movemap_cen_all()->get_bb(res) ) {
					pose.set_phi(res, -150.0);
					pose.set_psi(res, 150.0);
				}//if mm says ok
			}//for each residue
		}//if loop is extended
	}//for each loop

	return;
}


//option system replacement getters and setters
/// @brief run RMSD calculations
bool AnchoredDesignMover::get_rmsd() const {return rmsd_;}
/// @brief run only RMSD calculations against this input, don't do actual AnchoredDesign
std::string const & AnchoredDesignMover::get_RMSD_only_this() const {return RMSD_only_this_;}
/// @brief delete the input sidechains (independently from use_input_sc in the packer) - used to prevent leakage of sidechains in benchmarking mode
bool AnchoredDesignMover::get_delete_interface_native_sidechains() const {return delete_interface_native_sidechains_;}
/// @brief show_extended demonstrates that the code really forgets the input structure
bool AnchoredDesignMover::get_show_extended() const {return show_extended_;}
/// @brief randomize_input_sequence to complement loop extension in forgetting the input
bool AnchoredDesignMover::get_randomize_input_sequence() const {return randomize_input_sequence_;}
/// @brief pick a different cutpoint than the input; useful when you want to sample cutpoints
bool AnchoredDesignMover::get_vary_cutpoints() const {return vary_cutpoints_;}
/// @brief skip the perturbation step - useful when you already have a good structure
bool AnchoredDesignMover::get_refine_only() const {return refine_only_;}
/// @brief filter based on total complex score
core::Real AnchoredDesignMover::get_filter_score() const {return filter_score_;}
/// @brief filter based on complex SASA
core::Real AnchoredDesignMover::get_filter_SASA() const {return filter_SASA_;}
/// @brief filter based on omega angles in the loops - filter out cis omegas
bool AnchoredDesignMover::get_filter_omega() const {return use_filter_omega_;}
/// @brief whether to automatically initialize from the options system; defaults to true
bool AnchoredDesignMover::get_autoinitialize() const {return autoinitialize_;}

/// @brief run RMSD calculations
void AnchoredDesignMover::set_rmsd(bool const rmsd) { rmsd_ = rmsd;}
/// @brief run only RMSD calculations against this input, don't do actual AnchoredDesign
void AnchoredDesignMover::set_RMSD_only_this(std::string const & RMSD_only_this) { RMSD_only_this_ = RMSD_only_this;}
/// @brief delete the input sidechains (independently from use_input_sc in the packer) - used to prevent leakage of sidechains in benchmarking mode
void AnchoredDesignMover::set_delete_interface_native_sidechains(bool const delete_interface_native_sidechains) { delete_interface_native_sidechains_ = delete_interface_native_sidechains;}
/// @brief show_extended demonstrates that the code really forsets the input structure
void AnchoredDesignMover::set_show_extended(bool const show_extended) { show_extended_ = show_extended;}
/// @brief randomize_input_sequence to complement loop extension in forgetting the input
void AnchoredDesignMover::set_randomize_input_sequence(bool const randomize_input_sequence) { randomize_input_sequence_ = randomize_input_sequence;}
/// @brief pick a different cutpoint than the input { _ = ;} useful when you want to sample cutpoints
void AnchoredDesignMover::set_vary_cutpoints(bool const vary_cutpoints) { vary_cutpoints_ = vary_cutpoints;}
/// @brief skip the perturbation step - useful when you already have a good structure
void AnchoredDesignMover::set_refine_only(bool const refine_only) { refine_only_ = refine_only;}
/// @brief filter based on total complex score
void AnchoredDesignMover::set_filter_score(core::Real const filter_score) {
	filter_score_ = filter_score;
	use_filter_score_ = true;}
/// @brief filter based on complex SASA
void AnchoredDesignMover::set_filter_SASA(core::Real const filter_SASA) {
	filter_SASA_ = filter_SASA;
	use_filter_SASA_ = true;}
/// @brief filter based on omega angles in the loops - filter out cis omegas
void AnchoredDesignMover::set_filter_omega(bool const filter_omega) { use_filter_omega_ = filter_omega;}
/// @brief whether to automatically initialize from the options system; defaults to true
void AnchoredDesignMover::set_autoinitialize(bool const autoinitialize) { autoinitialize_ = autoinitialize;}

void AnchoredDesignMover::read_options(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys::AnchoredDesign;

	//bool rmsd_;
	rmsd_ = option[ OptionKeys::AnchoredDesign::rmsd ].value();

	//std::string RMSD_only_this_;
	if ( option[ testing::RMSD_only_this ].user() ) {
		RMSD_only_this_ = option[ testing::RMSD_only_this ].value();
	} else {
		RMSD_only_this_ = EMPTY_STRING;
	}

	//bool delete_interface_native_sidechains_;
	delete_interface_native_sidechains_ = option[ testing::delete_interface_native_sidechains].value();

	//bool show_extended_;
	show_extended_ = option[show_extended].value();

	//bool randomize_input_sequence_;
	randomize_input_sequence_ = option[ testing::randomize_input_sequence ].value();

	//bool vary_cutpoints_;
	vary_cutpoints_ = option[ vary_cutpoints ].value();

	//bool refine_only_;
	refine_only_ = option[ refine_only ].value();

	//core::Real filter_score_;
	if ( option[ OptionKeys::AnchoredDesign::filters::score].user() ) {
		filter_score_ = option[ OptionKeys::AnchoredDesign::filters::score].value();
		use_filter_score_ = true;
	} else {
		use_filter_score_ = false;
	}

	//core::Real filter_SASA_;
	if ( option[ OptionKeys::AnchoredDesign::filters::sasa].user() ) {
		filter_SASA_ = option[ OptionKeys::AnchoredDesign::filters::sasa].value();
		use_filter_SASA_ = true;
	} else {
		use_filter_SASA_ = false;
	}

	//bool filter_omega_;
	use_filter_omega_ = option[ OptionKeys::AnchoredDesign::filters::omega].value();

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////AnchoredPerturbMover////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
AnchoredPerturbMover::AnchoredPerturbMover( protocols::anchored_design::AnchorMoversDataOP interface_in) :
	Mover(), interface_(std::move( interface_in ))
{
	protocols::moves::Mover::type( "AnchoredPerturb" );
	read_options();
}

AnchoredPerturbMover::~AnchoredPerturbMover() = default;

/// @details AnchoredPerturbMover takes a pose, swaps it into centroid mode, and perturbs the structure of
///its mobile loops.  The perturbation step is under Monte Carlo control so it won't finish with a terrible
///structure.  It adds sidechains back on (unchanged relative to the alpha carbon) at the end, but does not
///repack those sidechains.
void AnchoredPerturbMover::apply( core::pose::Pose & pose )
{
	clock_t starttime = clock();
	T_perturb << "entering perturb steps" << std::endl;

	core::pose::Pose const saved_input_pose( pose ); //used to return sidechains later

	core::pose::PoseOP posecopy( nullptr );
	//if(debug_) posecopy = new core::pose::Pose( pose ); //after centroid-izing
	int counter(1);
	std::stringstream outputfilename;

	//score the protein
	T_perturb << "fullatom score of starting PDB: " << (*(interface_->get_fullatom_scorefunction()))(pose) << std::endl;
	(*(interface_->get_fullatom_scorefunction())).show( T_perturb, pose );
	T_perturb << std::flush; //show doesn't flush the buffer

	protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap(core::chemical::CENTROID);
	typeset_swap.apply( pose );
	if ( debug_ ) posecopy = core::pose::PoseOP( new core::pose::Pose( pose ) );

	//centroid score
	T_perturb << "centroid score of starting PDB: " << (*(interface_->get_centroid_scorefunction()))(pose) << std::endl;
	(*(interface_->get_centroid_scorefunction())).show( T_perturb, pose );
	T_perturb << std::flush; //show doesn't flush the buffer

	/*
	for each loop
	perturbmover of either type
	ccd close
	minimize
	MC evaluate
	*/

	//make the sequence mover.  other movers are inserted into it (their handles go out of scope!)
	protocols::moves::SequenceMoverOP perturb_sequence( new protocols::moves::SequenceMover() );
	protocols::moves::SequenceMoverOP allloops_subsequence( new protocols::moves::SequenceMover() );
	protocols::moves::RandomMoverOP backbone_mover( new protocols::moves::RandomMover() );

	//loop to add perturbing movers for each loop (this setup does each loop equally and sequentially)
	core::Size num_loops(interface_->num_loops());
	for ( core::Size i(1); i <= num_loops; ++i ) {
		///////////////////////////generate perturb mover///////////////////////////////
		//somewhat complex logic here, as there are three options for internal loops:
		// 1) SmallMover plus CCD closure
		// 2) Fragment insertion plus CCD closure, which is chosen over 1) if fragments are present
		// 3) Kinematic loop closure, which does not need separate perturb/close steps
		// We could also allow 3 with 1 OR 2 - try kinematic sometimes, and CCD sometimes
		// there are two possibilities for terminal flexible regions: SmallMover or fragments (no loop closure needed).

		protocols::moves::SequenceMoverOP oneloop_subsequence( new protocols::moves::SequenceMover() );
		protocols::moves::MoverOP perturb_mover;
		protocols::moves::RandomMoverOP oneloop_random( new protocols::moves::RandomMover() );
		//handles for clarity
		core::Size const loop_start(interface_->loop(i).start());
		core::Size const loop_end(interface_->loop(i).stop());
		bool const internal(!( pose.residue(loop_start).is_terminus() || pose.residue(loop_end).is_terminus() ));

		if ( perturb_CCD_off_ && perturb_KIC_off_ ) utility_exit_with_message( "cannot pass both AnchoredDesign::perturb_CCD_off AND AnchoredDesign::perturb_KIC_off; this turns off both types of loop remodeling" );

		//option 3: kinematic; no perturbation or closure needed
		if ( internal && !perturb_KIC_off_ ) {
			//make kinematic mover
			using protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
			KinematicMoverOP kin_mover( new KinematicMover() );//temperature, default 0.8
			KinematicMoverCAP kin_mover_cap(kin_mover);
			kin_mover->set_temperature( perturb_temp_ );

			//set up kinematic mover - this is borrowed from loops_main.cc, perturb_one_loop_with_alc(), SVN 24219, #932-946
			kin_mover->set_vary_bondangles( false );
			kin_mover->set_sample_nonpivot_torsions( nonpivot_torsion_sampling_ );
			kin_mover->set_rama_check( true );
			kin_mover->set_idealize_loop_first( false );

			//make kinematic perturber
			using protocols::loops::loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber;
			using protocols::loops::loop_closure::kinematic_closure::KinematicPerturberOP;
			KinematicPerturberOP TsamplingKP( new TorsionSamplingKinematicPerturber(kin_mover_cap) );
			TsamplingKP->set_movemap(interface_->movemap_cen(i));
			kin_mover->set_perturber(TsamplingKP);

			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapperOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapper;
			KinematicWrapperOP kin_wrapper( new KinematicWrapper(kin_mover, interface_->loop(i)) );
			kin_wrapper->respect_this_movemap(interface_->movemap_cen(i));
			oneloop_random->add_mover(kin_wrapper, 1);
			T_perturb << "creating KinematicWrapper with loop " << loop_start << " " << loop_end << std::endl;
		}//if using kinematic

		//if this is a terminal extension, or if we want CCD
		if ( !internal || !perturb_CCD_off_ ) {
			//now, check whether perturb is fragment based or not
			if ( !interface_->get_frags() /*if NULL*/ || no_frags_ ) {
				core::Size nmoves(5);
				using protocols::simple_moves::SmallMover;
				using protocols::simple_moves::SmallMoverOP;
				protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover(interface_->movemap_cen(i), 0.8, nmoves) );
				small_mover->angle_max( 'H', 180.0 );
				small_mover->angle_max( 'E', 180.0 );
				small_mover->angle_max( 'L', 180.0 );
				perturb_mover = small_mover;
				T_perturb << "creating SmallMover-based perturbation for loop " << loop_start << " " << loop_end << std::endl;
			} else { //smallmover based perturb
				/*T_perturb << "For loop start " << loop_start << " end " << loop_end
				<< " fragments exist; preferring fragments over small moves" << std::endl;*/
				using protocols::simple_moves::ClassicFragmentMover;
				protocols::simple_moves::ClassicFragmentMoverOP frag_mover( new ClassicFragmentMover(interface_->get_frags(), interface_->movemap_cen(i)) );
				frag_mover->enable_end_bias_check(false);
				perturb_mover = frag_mover;
				T_perturb << "creating fragment-based perturbation for loop " << loop_start << " " << loop_end << std::endl;
			} //fragment based perturb

			oneloop_subsequence->add_mover(perturb_mover);

			if ( debug_ ) {
				outputfilename.str(""); //clears stream
				outputfilename << "perturbed_" << counter;
				debug_dump_pose(pose, posecopy, outputfilename.str(), perturb_mover);
			}

			//now we need to add on closure, if the loop is not terminal and we didn't use KinematicMover above
			if ( internal && !perturb_CCD_off_ ) {
				///////////////////////////generate CCD close mover///////////////////////////////
				using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
				oneloop_subsequence->add_mover(moves::MoverOP( new CCDLoopClosureMover(interface_->loop(i), interface_->movemap_cen_omegafixed(i)) ));
				T_perturb << "creating CCD-closure after perturbation for loop " << loop_start << " " << loop_end << std::endl;

				if ( debug_ ) {
					outputfilename.str(""); //clears stream
					outputfilename << "perturbed_CCD_" << counter;
					debug_dump_pose(pose, posecopy, outputfilename.str(), oneloop_subsequence);
				}//debug
			}//if internal (for CCD)
			oneloop_random->add_mover(oneloop_subsequence, 1);
		}//if (CCD) or terminal

		++counter;

		allloops_subsequence->add_mover(oneloop_random); //keep track of all loops as a sequence
		backbone_mover->add_mover(oneloop_random, 1.0); //add to random mover with large weight
	} //end for loop over protein loops

	backbone_mover->add_mover(allloops_subsequence, 0.25); //run all loops combinatorially some small fraction of time
	//backbone_mover now randomly chooses one loop (equal probabilities) or all at once a small amount of the time
	perturb_sequence->add_mover(backbone_mover);

	if ( debug_ ) debug_dump_pose(pose, posecopy, "perturbed_nomin", allloops_subsequence);

	/////////////////////////minimizer mover/////////////////////////////////////////
	using protocols::simple_moves::MinMoverOP;
	using protocols::simple_moves::MinMover;
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(
		interface_->movemap_cen_all(),
		interface_->get_centroid_scorefunction_min(),
		min_type_,
		0.01,
		true /*use_nblist*/ ) );

	perturb_sequence->add_mover(min_mover);

	if ( debug_ ) {
		debug_dump_pose(pose, posecopy, "perturbed_min", min_mover);
		dump_cutpoint_info( pose );
		debug_dump_pose(pose, posecopy, "perturbed_all", perturb_sequence);
	}

	/////////////////////////wrap the sequence in a trial/////////////////////////////////////////
	//make the monte carlo object
	using protocols::moves::MonteCarloOP;
	using protocols::moves::MonteCarlo;
	MonteCarloOP mc( new MonteCarlo(
		pose,
		*(interface_->get_centroid_scorefunction()),
		perturb_temp_ ) );//temperature, default 0.8
	protocols::moves::TrialMoverOP trial_mover( new protocols::moves::TrialMover( perturb_sequence, mc ) );

	/////////////////////////wrap in a for loop output preference/////////////////////////////////
	T_perturb << "   Current     Low    total cycles =" << perturb_cycles_ << std::endl;
	for ( core::Size i = 1; i <= perturb_cycles_; ++i ) {
		trial_mover->apply( pose );
		T_perturb << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}//end the exciting for loop
	mc -> recover_low( pose );

	if ( perturb_show_ ) {
		using namespace protocols::jd2;
		JobDistributor::get_instance()->job_outputter()->other_pose(JobDistributor::get_instance()->current_job(), pose, "perturbed_centroid_final");
	}

	//show centroid score (duplicates last line above)
	T_perturb << "centroid score of final perturbed PDB: " << (*(interface_->get_centroid_scorefunction()))(pose)
		<< std::endl;
	(*(interface_->get_centroid_scorefunction())).show( T_perturb, pose );
	T_perturb << std::flush; //show doesn't flush the buffer

	if ( debug_ ) {
		using namespace protocols::jd2;
		JobDistributor::get_instance()->job_outputter()->other_pose(JobDistributor::get_instance()->current_job(), pose, "perturbed_prerefullatom");
	}
	protocols::simple_moves::ReturnSidechainMover return_sidechains( saved_input_pose );
	return_sidechains.apply( pose );

	/////////////////////////////generate full repack&minimize mover//////////////////////////////
	using core::pack::task::TaskFactoryOP; using core::pack::task::TaskFactory;
	TaskFactoryOP task_factory( new TaskFactory(*(interface_->get_late_factory())) ); //late factory = more rotamers
	task_factory->push_back(core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ));

	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
	pack_mover->task_factory( task_factory );
	pack_mover->score_function( interface_->get_fullatom_scorefunction() );

	//using protocols::simple_moves::MinMoverOP; //using protocols::simple_moves::MinMover;
	protocols::simple_moves::MinMoverOP min_mover_fa( new protocols::simple_moves::MinMover(
		interface_->movemap_cen_all(), //even though this is fullatom; we do not yet want to minimize inside the anchor if we are using constraints
		interface_->get_fullatom_scorefunction(),
		min_type_,
		0.01,
		true /*use_nblist*/ ) );

	//definitely want sidechain minimization here
	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover_fa( new protocols::simple_moves::TaskAwareMinMover(min_mover_fa, task_factory) );
	pack_mover->apply( pose );
	TAmin_mover_fa->apply( pose );

	//score the protein
	T_perturb << "fullatom score of perturbed, refullatomized, repacked/minimized PDB: "
		<< (*(interface_->get_fullatom_scorefunction()))(pose) << std::endl;
	(*(interface_->get_fullatom_scorefunction())).show( T_perturb, pose );
	T_perturb << std::flush; //show doesn't flush the buffer

	clock_t stoptime = clock();
	T_perturb << "One perturb took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;
	T_perturb << "perturb steps complete" << std::endl;

}//AnchoredPerturbMover::apply


std::string
AnchoredPerturbMover::get_name() const {
	return "AnchoredPerturbMover";
}

//option system replacement
/// @brief debugging mode activates a bunch of extra output
bool AnchoredPerturbMover::get_debug() const { return debug_;}
/// @brief do not perform CCD style closure (use KIC only)
bool AnchoredPerturbMover::get_perturb_CCD_off() const { return perturb_CCD_off_;}
/// @brief do not perform KIC style closure (use CCD only)
bool AnchoredPerturbMover::get_perturb_KIC_off() const { return perturb_KIC_off_;}
/// @brief use nonpivot torsion sampling for KIC?
bool AnchoredPerturbMover::get_nonpivot_torsion_sampling() const { return nonpivot_torsion_sampling_;}
/// @brief MC temperature
core::Real AnchoredPerturbMover::get_perturb_temp() const { return perturb_temp_;}
/// @brief number of MC cycles
core::Size AnchoredPerturbMover::get_perturb_cycles() const { return perturb_cycles_;}
/// @brief do not use fragments?
bool AnchoredPerturbMover::get_no_frags() const { return no_frags_;}
/// @brief what minimizer type to use?
std::string const & AnchoredPerturbMover::get_min_type() const { return min_type_;}
/// @brief show perturb result structure?
bool AnchoredPerturbMover::get_perturb_show() const { return perturb_show_;}

/// @brief debugging mode activates a bunch of extra output
void AnchoredPerturbMover::set_debug(bool const debug) { debug_= debug;}
/// @brief do not perform CCD style closure (use KIC only)
void AnchoredPerturbMover::set_perturb_CCD_off(bool const perturb_CCD_off) { perturb_CCD_off_= perturb_CCD_off;}
/// @brief do not perform KIC style closure (use CCD only)
void AnchoredPerturbMover::set_perturb_KIC_off(bool const perturb_KIC_off) { perturb_KIC_off_= perturb_KIC_off;}
/// @brief use nonpivot torsion sampling for KIC?
void AnchoredPerturbMover::set_nonpivot_torsion_sampling(bool const nonpivot_torsion_sampling) { nonpivot_torsion_sampling_= nonpivot_torsion_sampling;}
/// @brief MC temperature
void AnchoredPerturbMover::set_perturb_temp(core::Real const perturb_temp) { perturb_temp_= perturb_temp;}
/// @brief number of MC cycles
void AnchoredPerturbMover::set_perturb_cycles(core::Size const perturb_cycles) { perturb_cycles_= perturb_cycles;}
/// @brief do not use fragments?
void AnchoredPerturbMover::set_no_frags(bool const no_frags) { no_frags_= no_frags;}
/// @brief what minimizer type to use?
void AnchoredPerturbMover::set_min_type(std::string const & min_type) { min_type_= min_type;}
/// @brief show perturb result structure?
void AnchoredPerturbMover::set_perturb_show(bool const perturb_show) { perturb_show_= perturb_show;}

void AnchoredPerturbMover::read_options(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys::AnchoredDesign;

	//bool debug_;
	debug_ = option[ debug ].value();

	//bool perturb_CCD_off_;
	perturb_CCD_off_ = option[ perturb_CCD_off ].value();

	//bool perturb_KIC_off_;
	perturb_KIC_off_ = option[ perturb_KIC_off ].value();

	/// @brief use nonpivot torsion sampling for KIC?
	nonpivot_torsion_sampling_ = option[ OptionKeys::loops::nonpivot_torsion_sampling ].value();

	//core::Real perturb_temp_;
	perturb_temp_ = option[ perturb_temp ].value();

	//core::Size perturb_cycles_;
	perturb_cycles_ = option[ perturb_cycles ].value();

	//bool no_frags_;
	no_frags_ = option[ no_frags ].value();

	//std::string min_type_;
	min_type_ = option[ OptionKeys::run::min_type ].value();

	//bool perturb_show_;
	perturb_show_ = option[ perturb_show ].value();

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////AnchoredRefineMover/////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
AnchoredRefineMover::AnchoredRefineMover( protocols::anchored_design::AnchorMoversDataOP interface_in) :
	Mover(), interface_(std::move( interface_in ))
{
	protocols::moves::Mover::type( "AnchoredRefine" );
	read_options();
}

AnchoredRefineMover::~AnchoredRefineMover() = default;


/// @details
void AnchoredRefineMover::apply( core::pose::Pose & pose )
{
	clock_t starttime = clock();
	T_refine << "entering refine steps" << std::endl;

	//variables used for debugging output
	core::pose::PoseOP posecopy( nullptr );
	if ( debug_ ) posecopy = core::pose::PoseOP( new core::pose::Pose(pose) );
	int counter(1);
	//std::stringstream outputfilename;

	//score the protein
	(*(interface_->get_fullatom_scorefunction()))(pose);
	T_refine << "fullatom score upon entering refine mode:" << std::endl;
	(*(interface_->get_fullatom_scorefunction())).show( T_refine, pose );
	T_refine << std::flush; //show doesn't flush the buffer

	//refinement steps:
	/*
	first, do a repack/minimize of the interface (handles bad rotamers coming in from perturb mode)
	for each loop
	small backbone movements (small and shear)
	CCD
	rotamer trials
	minimize
	monte carlo check
	every so often, do a repack/minimize instead of backbone perturbations
	*/

	///////////////////////////////////////////MC object////////////////////////////////
	//make the MC object ahead of time
	using protocols::moves::MonteCarloOP;
	using protocols::moves::MonteCarlo;
	MonteCarloOP mc( new MonteCarlo(
		pose,
		*(interface_->get_fullatom_scorefunction()),
		refine_temp_ ) ); //temperature, default 0.8

	/////////////////////////////generate full repack mover//////////////////////////////
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
	pack_mover->task_factory( interface_->get_task_factory() );
	pack_mover->score_function( interface_->get_fullatom_scorefunction() );

	//////////////////////////////////generate minimizer mover/////////////////////////
	using protocols::simple_moves::MinMoverOP;
	using protocols::simple_moves::MinMover;
	protocols::simple_moves::MinMoverOP packing_min_mover( new protocols::simple_moves::MinMover(
		interface_->movemap_fa_all(),
		interface_->get_fullatom_scorefunction(),
		min_type_,
		0.01,
		true /*use_nblist*/ ) );

	using protocols::simple_moves::TaskAwareMinMoverOP;
	using protocols::simple_moves::TaskAwareMinMover;
	protocols::simple_moves::TaskAwareMinMoverOP packing_TAmin_mover( new protocols::simple_moves::TaskAwareMinMover(packing_min_mover, interface_->get_task_factory()) );

	////////////////////////////////////create repacking sequence///////////////////////////
	using protocols::moves::SequenceMover;
	using protocols::moves::SequenceMoverOP;
	SequenceMoverOP repack_sequence( new SequenceMover(pack_mover, packing_TAmin_mover) );

	//////////////////////////create refinement sequence////////////////////////////////
	//make the refine sequence mover.  other movers are inserted into it (their handles go out of scope!)
	protocols::moves::SequenceMoverOP refine_sequence( new protocols::moves::SequenceMover );
	protocols::moves::SequenceMoverOP allloops_subsequence( new protocols::moves::SequenceMover() );
	protocols::moves::RandomMoverOP backbone_mover( new protocols::moves::RandomMover() );

	//loop to add perturbing movers for each loop (this setup does each loop equally and sequentially)
	core::Size num_loops(interface_->num_loops());
	for ( core::Size i(1); i <= num_loops; ++i ) {

		///////////////////////////generate backbone mover///////////////////////////////
		//as before, we might want CCD, ALC, or both
		//we might have terminal loops

		protocols::moves::SequenceMoverOP oneloop_subsequence( new protocols::moves::SequenceMover() );
		protocols::moves::MoverOP perturb_mover;
		protocols::moves::RandomMoverOP oneloop_random( new protocols::moves::RandomMover() );
		//handles for clarity
		core::Size const loop_start(interface_->loop(i).start());
		core::Size const loop_end(interface_->loop(i).stop());
		bool const internal(!( pose.residue(loop_start).is_terminus() || pose.residue(loop_end).is_terminus() ));
		if ( refine_CCD_off_ && refine_KIC_off_ ) utility_exit_with_message( "cannot pass both AnchoredDesign::refine_CCD_off AND AnchoredDesign::refine_KIC_off; this turns off both types of loop remodeling" );

		//option 3: kinematic; no perturbation or closure needed
		if ( internal && !refine_KIC_off_ ) {
			//make kinematic mover
			using protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicMoverCAP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
			KinematicMoverOP kin_mover( new KinematicMover() );//temperature, default 0.8
			KinematicMoverCAP kin_mover_cap(kin_mover);
			kin_mover->set_temperature( refine_temp_ );

			//set up kinematic mover - this is borrowed from loops_main.cc, perturb_one_loop_with_alc(), SVN 24219, #932-946
			kin_mover->set_vary_bondangles( false );
			kin_mover->set_sample_nonpivot_torsions( nonpivot_torsion_sampling_ );
			kin_mover->set_rama_check( true );
			kin_mover->set_idealize_loop_first( false );

			//make kinematic perturber
			using protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber;
			using protocols::loops::loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP;
			VicinitySamplingKinematicPerturberOP TsamplingKP( new VicinitySamplingKinematicPerturber(kin_mover_cap) );
			TsamplingKP->set_movemap(interface_->movemap_fa(i));
			//TsamplingKP->set_sample_vicinity( vicinity_sampling_ );
			TsamplingKP->set_degree_vicinity( vicinity_degree_ );
			kin_mover->set_perturber(TsamplingKP);

			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapperOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapper;
			KinematicWrapperOP kin_wrapper( new KinematicWrapper(kin_mover, interface_->loop(i)) );
			kin_wrapper->respect_this_movemap(interface_->movemap_fa(i));
			oneloop_random->add_mover(kin_wrapper, 1);
			T_perturb << "creating KinematicWrapper with loop " << loop_start << " " << loop_end << std::endl;
		}//if using kinematic

		if ( !internal || !refine_CCD_off_ ) { //termini and ccd loops

			core::Size nmoves(1);
			protocols::moves::MoverOP small_mover( new protocols::simple_moves::SmallMover(interface_->movemap_fa(i), 0.8, nmoves) );
			protocols::moves::MoverOP shear_mover( new protocols::simple_moves::ShearMover(interface_->movemap_fa(i), 0.8, nmoves) );
			T_refine << "creating Small & ShearMover for loop " << loop_start << " " << loop_end << std::endl;

			//generate "real" movers
			//this will do either a shear or small move
			protocols::moves::RandomMoverOP smallshear_mover( new protocols::moves::RandomMover() );
			smallshear_mover->add_mover(small_mover);
			smallshear_mover->add_mover(shear_mover);

			oneloop_subsequence->add_mover(smallshear_mover);

			//if not a terminal loop
			if ( internal && !refine_CCD_off_ ) {
				using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
				using protocols::moves::MoverOP;
				MoverOP CCD_mover( new CCDLoopClosureMover(interface_->loop(i), interface_->movemap_fa_omegafixed(i)) );
				oneloop_subsequence->add_mover(CCD_mover);
				T_refine << "creating CCD-closure after perturbation for loop " << loop_start << " " << loop_end << std::endl;

				/*if(debug_){
				outputfilename.str(""); //clears stream
				outputfilename << "refine_smallmoved" << counter << ".pdb";
				small_mover->apply(*posecopy);
				posecopy->dump_pdb(outputfilename.str());

				outputfilename.str(""); //clears stream
				outputfilename << "refine_smallmoved_plus_CCD" << counter << ".pdb";
				debug_dump_pose(pose, posecopy, outputfilename.str(), CCD_mover);

				outputfilename.str(""); //clears stream
				outputfilename << "refine_shearmoved" << counter << ".pdb";
				shear_mover->apply(*posecopy);
				posecopy->dump_pdb(outputfilename.str());

				outputfilename.str(""); //clears stream
				outputfilename << "refine_shearmoved_plus_CCD" << counter << ".pdb";
				debug_dump_pose(pose, posecopy, outputfilename.str(), CCD_mover);
				}*/
			}//if needing CCD closure
			oneloop_random->add_mover(oneloop_subsequence, 1);
		}//if terminal or needing CCD closure
		++counter;
		if ( debug_ ) *posecopy = pose;

		allloops_subsequence->add_mover(oneloop_random); //keep track of all loops as a sequence
		backbone_mover->add_mover(oneloop_random, 1.0); //add to random mover with large weight
	} //end for loop over protein loops

	backbone_mover->add_mover(allloops_subsequence, 0.25); //run all loops combinatorially some small fraction of time
	//backbone_mover now randomly chooses one loop (equal probabilities) or all at once a small amount of the time
	refine_sequence->add_mover(backbone_mover);

	////////////////////rotamer trials//////////////////////
	//rotamer trials work only for repacking, no design
	using core::pack::task::TaskFactoryOP;
	using core::pack::task::TaskFactory;
	TaskFactoryOP rt_task_factory( new TaskFactory(*(interface_->get_task_factory())) ); //local copy so we can modify it
	rt_task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );

	using protocols::simple_moves::RotamerTrialsMoverOP;
	using protocols::simple_moves::EnergyCutRotamerTrialsMover;
	protocols::simple_moves::RotamerTrialsMoverOP rt_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover(
		interface_->get_fullatom_scorefunction(),
		rt_task_factory,
		mc,
		0.01 /*energycut*/ ) );

	if ( debug_ ) debug_dump_pose(pose, posecopy, "refine_rt.pdb", rt_mover);
	refine_sequence->add_mover(rt_mover);

	/////////////////////////minimizer mover/////////////////////////////////////////
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(
		interface_->movemap_fa_all(),
		interface_->get_fullatom_scorefunction(),
		min_type_,
		0.01,
		true /*use_nblist*/ ) );
	refine_sequence->add_mover( min_mover );

	if ( debug_ ) {
		debug_dump_pose(pose, posecopy, "refine_min.pdb", ( packing_TAmin_mover ));
		debug_dump_pose(pose, posecopy, "refine_all.pdb", refine_sequence);
		dump_cutpoint_info(pose);
	}

	//////////////////cycle mover controls when backbone perturbations occur versus repacking///////////////
	//this mover should do backbone refinement first <many> times it is called, then repack, then repeat
	protocols::moves::CycleMoverOP refine_repack_cycle( new protocols::moves::CycleMover );
	for ( core::Size i = 1; i < refine_repack_cycles_; i++ ) refine_repack_cycle->add_mover(refine_sequence);
	refine_repack_cycle->add_mover(repack_sequence);
	//could swap so that repack_trial is first and remove its explicit apply() above

	////////////////////////////////////Trial mover wraps sequence/////////////////////////////////////////
	protocols::moves::TrialMoverOP refine_master( new protocols::moves::TrialMover( refine_repack_cycle, mc ) );

	/////////////////////////wrap in a for loop output preference/////////////////////////////////

	core::Size const refine_latefactory = (refine_cycles_ -( refine_cycles_/3)); //magic number: final one third of cycles
	T_refine << "   Current     Low    total cycles =" << refine_cycles_
		<< ", second resfile at " << refine_latefactory << std::endl;
	for ( core::Size i = 1; i <= refine_cycles_; ++i ) {

		if ( i == refine_latefactory ) {
			pack_mover->task_factory(interface_->get_late_factory());
			//RT still needs restrict to repacking
			TaskFactoryOP rt_late_factory( new TaskFactory(*(interface_->get_late_factory())) );
			rt_late_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
			rt_mover->task_factory(rt_late_factory);
		}

		if ( i == refine_cycles_ ) {
			repack_sequence->apply(pose);
			mc->boltzmann(pose);
		} else {
			refine_master->apply( pose );
		}
		T_refine << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}//end the exciting for loop
	mc->recover_low( pose );

	//score the protein
	(*(interface_->get_fullatom_scorefunction()))(pose);
	T_refine << "fullatom score exiting refine mode:" << std::endl;
	(*(interface_->get_fullatom_scorefunction())).show( T_refine, pose );
	T_perturb << std::flush; //show doesn't flush the buffer

	clock_t stoptime = clock();
	T_refine << "One refine took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;
	T_refine << "refine steps complete" << std::endl;
}//AnchoredRefineMover::apply()

std::string
AnchoredRefineMover::get_name() const {
	return "AnchoredRefineMover";
}

//option system replacement
/// @brief debugging mode activates a bunch of extra output
bool AnchoredRefineMover::get_debug() const { return debug_;}
/// @brief do not perform CCD style closure (use KIC only)
bool AnchoredRefineMover::get_refine_CCD_off() const { return refine_CCD_off_;}
/// @brief do not perform KIC style closure (use CCD only)
bool AnchoredRefineMover::get_refine_KIC_off() const { return refine_KIC_off_;}
/// @brief use nonpivot torsion sampling for KIC?
bool AnchoredRefineMover::get_nonpivot_torsion_sampling() const { return nonpivot_torsion_sampling_;}
/// @brief KIC use vicinity sampling?
bool AnchoredRefineMover::get_vicinity_sampling() const { return vicinity_sampling_;}
/// @brief KIC vicinity sampling degrees
core::Real AnchoredRefineMover::get_vicinity_degree() const { return vicinity_degree_;}
/// @brief MC temperature
core::Real AnchoredRefineMover::get_refine_temp() const { return refine_temp_;}
/// @brief number of MC cycles
core::Size AnchoredRefineMover::get_refine_cycles() const { return refine_cycles_;}
/// @brief what minimizer type to use?
std::string const & AnchoredRefineMover::get_min_type() const { return min_type_;}
/// @brief how many cycles between repack/design opportunities?
core::Size AnchoredRefineMover::get_refine_repack_cycles() const { return refine_repack_cycles_;}

/// @brief debugging mode activates a bunch of extra output
void AnchoredRefineMover::set_debug(bool const debug) { debug_= debug;}
/// @brief do not perform CCD style closure (use KIC only)
void AnchoredRefineMover::set_refine_CCD_off(bool const refine_CCD_off) { refine_CCD_off_= refine_CCD_off;}
/// @brief do not perform KIC style closure (use CCD only)
void AnchoredRefineMover::set_refine_KIC_off(bool const refine_KIC_off) { refine_KIC_off_= refine_KIC_off;}
/// @brief use nonpivot torsion sampling for KIC?
void AnchoredRefineMover::set_nonpivot_torsion_sampling(bool const nonpivot_torsion_sampling) { nonpivot_torsion_sampling_= nonpivot_torsion_sampling;}
/// @brief KIC use vicinity sampling?
void AnchoredRefineMover::set_vicinity_sampling(bool const vicinity_sampling) { vicinity_sampling_= vicinity_sampling;}
/// @brief KIC vicinity sampling degrees
void AnchoredRefineMover::set_vicinity_degree(core::Size const vicinity_degree) { vicinity_degree_= vicinity_degree;}
/// @brief MC temperature
void AnchoredRefineMover::set_refine_temp(core::Real const refine_temp) { refine_temp_= refine_temp;}
/// @brief number of MC cycles
void AnchoredRefineMover::set_refine_cycles(core::Size const refine_cycles) { refine_cycles_= refine_cycles;}
/// @brief what minimizer type to use?
void AnchoredRefineMover::set_min_type(std::string const & min_type) { min_type_= min_type;}
/// @brief how many cycles between repack/design opportunities?
void AnchoredRefineMover::set_refine_repack_cycles(core::Size const refine_repack_cycles) { refine_repack_cycles_= refine_repack_cycles;}

void AnchoredRefineMover::read_options(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys::AnchoredDesign;

	//bool debug_;
	debug_ = option[ debug ].value();

	//bool refine_CCD_off_;
	refine_CCD_off_ = option[ refine_CCD_off ].value();

	//bool refine_KIC_off_;
	refine_KIC_off_ = option[ refine_KIC_off ].value();

	/// bool nonpivot_torsion_sampling_;
	nonpivot_torsion_sampling_ = option[ OptionKeys::loops::nonpivot_torsion_sampling ].value();

	//bool vicinity_sampling_;
	vicinity_sampling_ = option[ OptionKeys::loops::vicinity_sampling ].value();

	//core::Real vicinity_degree_;
	vicinity_degree_ = option[ OptionKeys::loops::vicinity_degree ].value();

	//core::Real refine_temp_;
	refine_temp_ = option[ refine_temp ].value();

	//core::Size refine_cycles_;
	refine_cycles_ = option[ refine_cycles ].value();

	//std::string min_type_;
	min_type_ = option[ OptionKeys::run::min_type ].value();

	//core::Size refine_repack_cycles_;
	refine_repack_cycles_ = option[ refine_repack_cycles ].value();

	return;
}

}//AnchoredDesign
}//protocols
