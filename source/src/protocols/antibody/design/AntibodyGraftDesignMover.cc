// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyGraftDesignMover.cc
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDesignMoverGenerator.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>
#include <protocols/antibody/design/util.hh>

// Core Includes
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>

// Protocol Includes
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/util.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

// Numeric Includes
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

//Utility
#include <math.h>
#include <basic/Tracer.hh>
#include <utility/py/PyAssert.hh>
#include <map>
#include <boost/algorithm/string.hpp>


static thread_local basic::Tracer TR("antibody.design.AntibodyGraftDesignMover");

namespace protocols{
namespace antibody{
namespace design{
	using namespace protocols::antibody;
	using namespace protocols::grafting;
	using namespace protocols::antibody::clusters;
	using namespace protocols::antibody::constraints;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using core::Size;

AntibodyGraftDesignMover::AntibodyGraftDesignMover(AntibodyInfoOP ab_info):
	ab_info_(ab_info),
	graft_mover_(/* NULL */),
	scorefxn_(/* NULL */)
{
	overhang_ = 3;
	modeler_ = AntibodyDesignModelerOP( new AntibodyDesignModeler(ab_info_) );
	mover_generator_ = AntibodyDesignMoverGeneratorOP( new AntibodyDesignMoverGenerator(ab_info_) );

	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}
	read_command_line_options();
	set_defaults();
}

AntibodyGraftDesignMover::~AntibodyGraftDesignMover(){}

void
AntibodyGraftDesignMover::set_defaults(){
	///Conservative defaults.  Defaults here are also read and set from a database file.  To allow run-time manipulations and testing.
	print_tracer_info_ = true;
	set_post_graft_modeling_cycles(1);
	paratope_cdrs_.clear();
	paratope_cdrs_.resize(6, true);
	epitope_residues_.clear();
	cdr_set_options_.clear();
	cdr_graft_design_options_.clear();
	cdr_seq_design_options_.clear();
}

void
AntibodyGraftDesignMover::setup_options_classes(){
	if (cdr_set_options_.size() == 0){
		cdr_set_options_  = protocols::antibody::design::get_cdr_set_options();
	}
	if (cdr_graft_design_options_.size() == 0){
		cdr_graft_design_options_ = protocols::antibody::design::get_graft_design_options();
	}
	if (cdr_seq_design_options_.size() == 0){
		cdr_seq_design_options_ = protocols::antibody::design::get_seq_design_options();
	}

	///Set it up so that if we design the cdr, we load CDRs from the database.
	for (core::Size i = 1; i <= 6; ++i){
		CDRGraftDesignOptionsOP options = cdr_graft_design_options_[i];
		cdr_set_options_[i]->load(options->design());
	}

}

void
AntibodyGraftDesignMover::read_command_line_options(){
	set_graft_rounds(basic::options::option [basic::options::OptionKeys::antibody::design::graft_rounds]());
	set_keep_top_designs(basic::options::option [basic::options::OptionKeys::antibody::design::top_graft_designs]());
	set_dock_post_graft(basic::options::option [basic::options::OptionKeys::antibody::design::post_graft_dock]());
	set_rb_min_post_graft(basic::options::option [basic::options::OptionKeys::antibody::design::post_graft_rb_min]());
	set_dock_rounds(basic::options::option [basic::options::OptionKeys::antibody::design::dock_rounds]());
	initial_perturb_ = basic::options::option [basic::options::OptionKeys::antibody::design::initial_perturb] ();
	use_deterministic_algorithm_ = basic::options::option [basic::options::OptionKeys::antibody::design::use_deterministic]();
	benchmark_ = basic::options::option [basic::options::OptionKeys::antibody::design::benchmark_graft_design]();
	use_light_chain_type_ = basic::options::option [basic::options::OptionKeys::antibody::design::use_light_chain_type]();
	use_epitope_constraints_ = basic::options::option [basic::options::OptionKeys::antibody::design::use_epitope_constraints]();
	adapt_graft_ = basic::options::option [basic::options::OptionKeys::antibody::design::adapt_graft]();
	design_post_graft_ = basic::options::option [basic::options::OptionKeys::antibody::design::post_graft_seq_design]();

}

void
AntibodyGraftDesignMover::set_scorefunction(ScoreFunctionOP scorefxn){
	scorefxn_=scorefxn->clone();
}

void
AntibodyGraftDesignMover::set_keep_top_designs(core::Size top_designs){
	num_top_designs_ = top_designs;
}

void
AntibodyGraftDesignMover::set_dock_post_graft(bool dock_post_graft){
	dock_post_graft_ = dock_post_graft;
}

void
AntibodyGraftDesignMover::set_dock_rounds(core::Size dock_rounds){
	dock_rounds_ = dock_rounds;
}

void
AntibodyGraftDesignMover::set_rb_min_post_graft(bool rb_min_post_graft){
	rb_min_post_graft_ = rb_min_post_graft;
}
void
AntibodyGraftDesignMover::set_graft_rounds(core::Size graft_rounds){
	graft_rounds_=graft_rounds;
}

void
AntibodyGraftDesignMover::set_post_graft_modeling_cycles(core::Size cycles) {
	post_graft_modeling_cycles_ = cycles;
}


void
AntibodyGraftDesignMover::set_paratope_cdrs(const vector1<bool>& cdrs) {
	if (cdrs.size() != 6){
		utility_exit_with_message("Cannot setup paratope cdrs - number passed does not equal total number of cdrs!");
	}
	paratope_cdrs_ = cdrs;

}

void
AntibodyGraftDesignMover::set_epitope_residues(vector1<PDBNumbering> epitope_residues) {
	epitope_residues_ = epitope_residues;
}



void
AntibodyGraftDesignMover::setup_native_clusters(core::pose::Pose & pose){
	ab_info_->setup_CDR_clusters(pose);
}

void
AntibodyGraftDesignMover::set_cdr_set(CDRSet& cdr_set, core::Size overhang){
	cdr_set_ = cdr_set;
	overhang_ = overhang;
}

void
AntibodyGraftDesignMover::initialize_cdr_set(core::pose::Pose const & pose){
	AntibodyDatabaseManager manager = AntibodyDatabaseManager(ab_info_);
	cdr_set_ = manager.load_cdr_poses(cdr_set_options_, pose, use_light_chain_type_);
}

void
AntibodyGraftDesignMover::set_cdr_set_options(AntibodyCDRSetOptions cdr_set_options){
	cdr_set_options_ = cdr_set_options;
}

void
AntibodyGraftDesignMover::set_graft_design_options(AntibodyCDRGraftDesignOptions graft_design_options){
	cdr_graft_design_options_ = graft_design_options;
}

void
AntibodyGraftDesignMover::set_seq_design_options(AntibodyCDRSeqDesignOptions seq_design_options){
	cdr_seq_design_options_ = seq_design_options;
}

void
AntibodyGraftDesignMover::setup_default_graft_settings(){

	//CCDEndsGraftMover
	graft_mover_->set_cycles(20);
	graft_mover_->set_scaffold_flexibility(2, 2);
	graft_mover_->set_insert_flexibility(0, 0);
	graft_mover_->final_repack(true);
	graft_mover_->stop_at_closure(true);
	graft_mover_->neighbor_dis(mover_generator_->neighbor_detection_dis());
	graft_mover_->set_fa_scorefunction(scorefxn_);

	//graft_mover_->set_skip_sampling(true);

	//AnchoredGraftMover
	anchored_graft_mover_->set_cycles(100);
	anchored_graft_mover_->set_scaffold_flexibility(2, 2);
	anchored_graft_mover_->set_insert_flexibility(1, 1);
	anchored_graft_mover_->final_repack(true);
	anchored_graft_mover_->stop_at_closure(true);
	//anchored_graft_mover_->conservative_smallmover(true);
	anchored_graft_mover_->neighbor_dis(mover_generator_->neighbor_detection_dis());
	anchored_graft_mover_->set_fa_scorefunction(scorefxn_);
}

void
AntibodyGraftDesignMover::setup_scorefxn() {

	if (!scorefxn_){
		scorefxn_ = core::scoring::get_score_function(true);
		scorefxn_->apply_patch_from_file("antibody_design");
	}
}

void
AntibodyGraftDesignMover::setup_paratope_epitope_constraints(core::pose::Pose & pose){

	if (use_epitope_constraints_){
		paratope_epitope_cst_mover_->set_defaults(); //Clear everything
		paratope_epitope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_); //Regenerates paratope residues
		paratope_epitope_cst_mover_->constrain_to_epitope_residues(epitope_residues_, pose);
		paratope_epitope_cst_mover_->set_interface_distance(mover_generator_->interface_detection_dis());
		paratope_epitope_cst_mover_->apply(pose);
	}
	else {
		paratope_cst_mover_->set_defaults(); //Clear everything
		paratope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_);
		paratope_cst_mover_->set_interface_distance(mover_generator_->interface_detection_dis());
		paratope_cst_mover_->apply(pose);
	}
}

void
AntibodyGraftDesignMover::setup_random_start_pose(pose::Pose& pose, vector1<CDRNameEnum>& cdrs_to_design){

	TR << "Randomizing starting CDRs" << std::endl;
	print_tracer_info_ = false;
	core::Size max_attempts = 1000;
	bool original_adaptation = adapt_graft_;
	for (core::Size i = 1; i <= cdrs_to_design.size(); ++i){
		bool graft_closed = false;
		core::Size attempts = 0;
		CDRNameEnum cdr = cdrs_to_design[i];
		while (! graft_closed){
			core::Size cdr_index = numeric::random::rg().random_range(1, cdr_set_[cdr].size());
			graft_closed = this->apply_to_cdr(pose, cdr, cdr_index, false); //No minimization
			++attempts;

			//Infinite loop prevention for stubborn CDRs
			if (attempts == max_attempts){
				if (! adapt_graft_){
					adapt_graft_ = true;
					attempts = 0;
				}
				else{
					utility_exit_with_message(
						"Could not setup randomized starting pose as graft for " + ab_info_->get_CDR_name(cdr)+" Max attempts reached!");
				}
			}
		}
		TR << "Attempts "<<attempts << " " << ab_info_->get_CDR_name(cdr) << std::endl;
		adapt_graft_ = original_adaptation;
	}
	protocols::jd2::JobOP job = protocols::jd2::JobDistributor::get_instance()->current_job();
	std::string prefix = "initial_benchmark_perturbation";
	protocols::jd2::JobDistributor::get_instance()->job_outputter()->other_pose(job, pose, prefix);
	print_tracer_info_ = true;
}


bool
AntibodyGraftDesignMover::apply_to_cdr(pose::Pose & pose, CDRNameEnum cdr, core::Size index, bool min_post_graft /*true*/){
	using namespace protocols::docking;

	//Graft CDR into scaffold pose.  If unsuccessful, do not continue.
	pose.remove_constraints(); //Remove constraints to increase speed of insertion. Quicker to remove and re-add.
	CDRPose & current_cdr_pose = cdr_set_[cdr][index];

	if (! graft_in_cdr(pose, cdr, current_cdr_pose)) return false;
	setup_paratope_epitope_constraints(pose);

	//Setup the clusters in AntibodyInfo before any minimization can presumably change the cluster.
	//This is primarily for SeqDesign during minimization.
	//CDRClusterSet should presumably be part of the Pose DataCache, however, we will need AntibodyInfo's clone methods.
	//We will also presumably, want a more basic CDRClusterSet instead, which would basically be a vector of CDRClusterOps.

	core::Size grafted_cdr_length = ab_info_->get_CDR_length(cdr, pose);

	CDRClusterOP cluster( new CDRCluster(
			pose,
			cdr,
			grafted_cdr_length,
			current_cdr_pose.cluster,
			ab_info_->get_CDR_start(cdr, pose),
			current_cdr_pose.distance) );

	//May be a bit expensive copying wise, but very useful nonetheless.
	ab_info_->set_CDR_cluster(cdr, cluster);
	ab_info_->get_CDR_cluster_set()->set_cacheable_cluster_data_to_pose(pose);

	//This probably helps now with the paratope/epitope constraints.  Should it be outside of the graft step though?
	if (initial_perturb_){
		TR <<"Perturbing structure" << std::endl;
		std::string dock_chains = ab_info_->get_antibody_chain_string()+"_A";
		core::kinematics::FoldTree ft = pose.fold_tree();
		vector1< int > movable_jumps(1, 1);
		setup_foldtree(pose, dock_chains, movable_jumps);
		DockingInitialPerturbationOP perturber( new DockingInitialPerturbation(1, true /* slide */) );

		perturber->apply(pose);
		pose.fold_tree(ft);
	}


	if (min_post_graft){
		TR << "Beginning post-graft min" << std::endl;
		protocols::moves::MonteCarlo mc = protocols::moves::MonteCarlo(pose, *scorefxn_, 2.0);
		for (core::Size i = 1; i <= post_graft_modeling_cycles_; ++i){

			TR << "post graft modeling round: "<< i << std::endl;
			run_post_graft_min(pose, mc, cdr);
		}
		mc.recover_low(pose);
	}
	TR << "Graft complete" << std::endl;
	return true;
}

bool
AntibodyGraftDesignMover::graft_in_cdr(pose::Pose& pose, CDRNameEnum const cdr, CDRPose & cdr_pose){

	TR << "Grafting CDR from cluster " << ab_info_->get_cluster_name(cdr_pose.cluster) << " fragment "<< cdr_pose.pdb << std::endl;

	core::Size start = ab_info_->get_CDR_start(cdr, pose)-1;
	core::Size end = ab_info_->get_CDR_end(cdr, pose)+1;
	graft_mover_->set_insert_region(start, end);
	anchored_graft_mover_->set_insert_region(start, end);
	if (cdr==h3){
		graft_mover_->set_cycles(40);
	}
	pose::Pose temp_pose = pose;
	std::pair<bool, core::Size> cb = run_graft(temp_pose, cdr, cdr_pose, graft_mover_);
	if (cb.first && adapt_graft_){
		TR << "Graft not closed. Adapting graft closure" << std::endl;
		pose::Pose temp_pose = pose;
		cb = run_graft(temp_pose, cdr, cdr_pose, anchored_graft_mover_);
	}
	if (! cb.first){
		pose = temp_pose;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Print data about graft closure to parse from output. Return if not closed.

	if (print_tracer_info_){
		CDRClusterEnum cluster = cdr_pose.cluster;
		//core::id::AtomID c_0 = core::id::AtomID(pose.residue(cb.second).atom_index("C"), cb.second);
		//core::id::AtomID n_1 = core::id::AtomID(pose.residue(cb.second).atom_index("N"), cb.second + 1);
		//core::Real c_n_dis =pose.conformation().bond_length(c_0, n_1);

		TR.Trace  << "DATA GRAFT_CLOSURE "<< std::boolalpha << !cb.first <<" " << pose.pdb_info()->name() << " "
			<< ab_info_->get_CDR_name(cdr)<<" "<< ab_info_->get_cluster_name(cluster) <<" "<< cdr_pose.pdb << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (cb.first){

		return false;
	}
	else{

		return true;
	}
}

std::pair<bool, core::Size>
AntibodyGraftDesignMover::run_graft(pose::Pose& pose, CDRNameEnum const cdr, CDRPose & cdr_pose, AnchoredGraftMoverOP grafter){
	grafter->set_piece(*cdr_pose.pose, overhang_, overhang_);
	grafter->apply(pose);
	pose.pdb_info()->copy(*(cdr_pose.pose->pdb_info()), 1 + overhang_, cdr_pose.pose->total_residue()-overhang_, grafter->start()+1);
	pose.pdb_info()->obsolete(false);
	protocols::antibody::constraints::add_harmonic_cluster_cst_or_dihedral_cst(ab_info_, pose, cdr, cdr_pose.cluster);
	scorefxn_->score(pose); //Segfault prevention.

	modeler_->set_cdr_only(cdr, true);
	mover_generator_->set_cdr_only(cdr, true);

	TR << "Checking Geometry" << std::endl;
	std::pair<bool, core::Size> cb = check_cb(pose, grafter->get_loops());

	if (grafter->get_name()=="AnchoredGraftMover"){
		//Anchored graft mover can severely screw up the graft.  Attempt fix with cartmin.
		//modeler_->minimize_cdrs(pose, true/*min_sc*/, true /*neighbor_sc*/, false /*min_interface*/, true /*cartmin*/, false /*talaris_cartmin*/);
		mover_generator_->set_as_mover(true);
		mover_generator_->set_task_factory(NULL); //Reset any TaskFactory.
		mover_generator_->set_min_interface(false);
		mover_generator_->set_cartmin(true);
		mover_generator_->set_min_sc(true);
		mover_generator_->set_include_neighbor_sc(true);
		mover_generator_->generate_minimizer(pose);
		mover_generator_->apply(pose);

		cb = check_cb(pose, grafter->get_loops());
	}
	else if (cb.first){
		TR <<"attempting cb fix with bb cartesian min" << std::endl;
		//modeler_->minimize_cdrs(pose, false /*min_sc*/, false /*neighbor_sc*/, false /*min_interface*/, true /*cartmin*/, false /*talaris_cartmin*/);
		mover_generator_->set_as_mover(true);
		mover_generator_->set_task_factory(NULL); //Reset any TaskFactory.
		mover_generator_->set_min_interface(false);
		mover_generator_->set_cartmin(true);
		mover_generator_->set_min_sc(false);
		mover_generator_->set_include_neighbor_sc(false);
		mover_generator_->generate_minimizer(pose);
		mover_generator_->apply(pose);

		cb = check_cb(pose, grafter->get_loops());
	}
	return cb;
}

void
AntibodyGraftDesignMover::run_post_graft_min(pose::Pose& pose, protocols::moves::MonteCarlo & mc, CDRNameEnum cdr) {


	//Setup Neighbor CDR minimization - instruction file/options controlled..
	CDRGraftDesignOptionsOP options = cdr_graft_design_options_[cdr];

	vector1<bool> cdrs_to_min(6, false);
	utility::vector1<CDRNameEnum> neighbor_min = options->neighbor_min();

	//Convert from vector of cdrs to vector1 bool.
	cdrs_to_min[cdr] = true;
	for (core::Size i = 1; i <= neighbor_min.size(); ++i){
		cdrs_to_min[neighbor_min[i]] = true;
	}

	modeler_->set_cdrs(cdrs_to_min);
	mover_generator_->set_cdrs(cdrs_to_min);

	mover_generator_->set_as_mover(true);
	mover_generator_->set_min_sc(options->min_sc());
	mover_generator_->set_include_neighbor_sc(options->min_neighbor_sc());
	mover_generator_->set_min_interface(false);
	//mover_generator_->set_min_interface(options->min_rb()); //Not yet work with CDR min due to FT!

	//Setup design task factory.
	TaskFactoryOP design_tf;
	if (design_post_graft_){
		//Need to setup correct clusters for all CDRs for correct sequences. - Should be reflected in AntibodyInfoCOP of seq_design_creator
		design_tf = seq_design_creator_->generate_tf_seq_design_graft_design(pose, cdr, cdrs_to_min);
		mover_generator_->set_task_factory(design_tf);

		//Debug:
		core::pack::task::PackerTaskOP task = design_tf->create_task_and_apply_taskoperations(pose);
		task->show();

	}
	//Dock - before min?
	if (dock_post_graft_){
		for (core::Size i = 1; i<= dock_rounds_; i++){
			//Control this from cmd line
			TR << "Dock round "<< i <<std::endl;
			if (design_post_graft_){
				TaskFactoryOP dock_design_tf = seq_design_creator_->generate_tf_seq_design(pose);
				seq_design_creator_->disable_design_for_non_designing_cdrs(dock_design_tf, pose);
				dock_design_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface(1, mover_generator_->interface_detection_dis()) ));
				mover_generator_->set_task_factory(dock_design_tf);

				//Debug:

			}

			mover_generator_->generate_dock_low_res(pose);
			mover_generator_->apply(pose);

			mover_generator_->generate_repack_antigen_ab_interface(pose);
			mover_generator_->apply(pose);

			modeler_->minimize_interface(pose); //Seems to help.  mover_generator can't yet set this up.

			mover_generator_->generate_repack_antigen_ab_interface(pose);
			mover_generator_->apply(pose);

			mover_generator_->generate_dock_high_res(pose);
			mover_generator_->apply(pose);
			mc.boltzmann(pose); //Low res dock will definitely screw up the pose - so do MC after both low and high res
		}
	}
	//Repack CDR?
	if (options->mintype() == repack){
		//modeler_->repack_cdrs(pose, options->min_neighbor_sc());
		mover_generator_->generate_repack_cdrs(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);
	}
	//Relax CDR?
	else if (options->mintype() == relax){
		//modeler_->relax_cdrs(pose, options->min_neighbor_sc(), false, options->min_rb());
		mover_generator_->set_dualspace(false);
		mover_generator_->generate_relax(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);
	}
	else if (options->mintype() == dualspace) {
		//modeler_->relax_cdrs(pose, options->min_neighbor_sc(), false, options->min_rb(), true /*dualspace*/);
		mover_generator_->set_dualspace(true);
		mover_generator_->generate_relax(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);
	}
	//Minimize CDR?
	else if (options->mintype() == minimize){

		mover_generator_->generate_repack_cdrs(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);

		mover_generator_->set_cartmin(false);
		mover_generator_->generate_minimizer(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);
	}
	else if (options->mintype() == minimize_cartesian) {

		mover_generator_->generate_repack_cdrs(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);

		mover_generator_->set_cartmin(true);
		mover_generator_->generate_minimizer(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose);
	}
	//Interface Rigid Body minimize?
	if (rb_min_post_graft_){
		modeler_->minimize_interface(pose, true /* min_interface_sc */);
		mc.boltzmann(pose);
	}

	//Dock - post-min?
	if (dock_post_graft_){
		modeler_->minimize_interface(pose, true); //Seems to help
		mover_generator_->generate_dock_high_res(pose);
		mover_generator_->apply(pose);
		mc.boltzmann(pose); //Low res dock will definitely screw up the pose - so do MC after both low and high res
	}
}

void
AntibodyGraftDesignMover::check_for_top_designs(pose::Pose & pose){

	//Can be refactored to use utility::TopScoreSelector
	//From mc algorithm, you can have multiple poses that are equivalent...
	core::Real score = (*scorefxn_)(pose);
	vector1<core::Real>::iterator score_it = top_scores_.begin();
	vector1<pose::PoseOP>::iterator pose_it = top_designs_.begin();
	if (top_scores_.size()==0){
		top_scores_.push_back(score);
		top_designs_.push_back(core::pose::PoseOP( new Pose() ));
		*(top_designs_[top_designs_.size()]) = pose;
	}
	else{
		bool inserted = false;
		for (core::Size i = 1; i<=top_scores_.size(); ++i){
			if (score <= top_scores_[i]){
				top_scores_.insert(score_it+i-1, score);
				top_designs_.insert(pose_it+i-1, core::pose::PoseOP( new Pose() ));
				*(top_designs_[i]) = pose;
				inserted = true;
				break;
			}
		}
		if (! inserted && top_scores_.size() < num_top_designs_){
			top_scores_.push_back(score);
			top_designs_.push_back(core::pose::PoseOP( new Pose() ));
			*(top_designs_[top_designs_.size()]) = pose;
		}
		else if ( inserted && top_scores_.size() > num_top_designs_){
			top_scores_.pop_back();
			top_designs_.pop_back();
		}
	}
	//mc_->eval_lowest_score_pose(pose, false, true);
}

///@brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
vector1< vector1 < Size > >
AntibodyGraftDesignMover::get_cdr_set_index_list(){
	vector1< vector1< Size > > index_list;
	vector1<core::Size> cdr_set_totals(6, 0);

	for (core::SSize i=CDRNameEnum_start; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		cdr_set_totals[i] = cdr_set_[cdr].size();
	}

	vector1< Size > dummy_index(6, 0);
	get_all_graft_permutations(cdr_set_totals, index_list, dummy_index, 1);
	return index_list;

}

void
AntibodyGraftDesignMover::run_basic_mc_algorithm(pose::Pose& pose, vector1<CDRNameEnum>& cdrs_to_design){

	TR << "Running basic monte carlo algorithm " << std::endl;

	top_scores_.push_back((*scorefxn_)(pose));
	top_designs_.push_back(core::pose::PoseOP( new Pose() ));
	*(top_designs_[1]) = pose;
	mc_->set_last_accepted_pose(pose);

	//Setup Weights
	utility::vector1<core::Real> cdr_weights;
	numeric::random::WeightedSampler sampler;
	for (core::Size i = 1; i <= cdrs_to_design.size(); ++i){
		core::Real weight = cdr_graft_design_options_[cdrs_to_design[ i ]]->weight();
		
		//CDRNameEnum cdr = static_cast<CDRNameEnum>(cdrs_to_design[ i ]);
		//TR << ab_info_->get_CDR_name(cdr) <<" Weight " << weight << std::endl;
		cdr_weights.push_back(weight);
	}
	sampler.weights(cdr_weights);

	//Choose random weighted CDR, graft in CDR, minimize
	for (core::Size i = 1; i <= graft_rounds_; ++i){
		TR << "Graft round: " << i <<std::endl;
		CDRNameEnum cdr_type = cdrs_to_design[sampler.random_sample(numeric::random::rg())];
		core::Size cdr_index = numeric::random::rg().random_range(1, cdr_set_[cdr_type].size());

		bool graft_successful = apply_to_cdr(pose, cdr_type, cdr_index);
		if ( graft_successful ){
			check_for_top_designs(pose);
			mc_->boltzmann(pose);
		}
		else{
			pose = mc_->last_accepted_pose();
		}
	}
	mc_->recover_low(pose);
	mc_->show_counters();
}

void
AntibodyGraftDesignMover::run_deterministic_graft_algorithm(pose::Pose & pose,vector1<CDRNameEnum>& cdrs_to_design){
	//Note:  Not feasible with >= 4 CDRs to try at each position. Should not be used with docking on.
	//Temporary fix to deterministically graft one CDR.  Needs to be fixed correctly soon.
	CDRNameEnum cdr = cdrs_to_design[1];
	for (core::Size i = 1; i <= cdr_set_[cdr].size(); ++i){
		TR << "Graft round: " << i << std::endl;
		pose::Pose trial_pose = pose;
		bool graft_successful = apply_to_cdr(trial_pose, cdr, i);
		if (! graft_successful){
			TR << "Graft unsuccessful.  Skipping combination..." << std::endl;
			continue;
		}
		//Only check for top designs after all combinations have been grafted.
		check_for_top_designs(trial_pose);
		mc_->eval_lowest_score_pose(trial_pose, false, true);
	}
	mc_->recover_low(pose);
	mc_->show_counters();
}

void
AntibodyGraftDesignMover::show(std::ostream & output) const{
	//Will go in SHOW.

	utility::vector1<std::string> strategies(SeqDesignStrategyEnum_total, "");
	strategies[ seq_design_profiles ] = "Profile-Based";
	strategies[ seq_design_conservative ] = "Conservative";
	strategies[ seq_design_basic ] = "Basic";

	output << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
	for (core::Size i=1; i<=6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		CDRSetOptionsCOP set_options = cdr_set_options_[cdr];
		CDRGraftDesignOptionsCOP options = cdr_graft_design_options_[cdr];
		CDRSeqDesignOptionsCOP seq_options = cdr_seq_design_options_[cdr];

		output << "//////// " << ab_info_->get_CDR_name(cdr) << " ///////////////////////////////////////////////"<< std::endl;
		output << "//// "<< std::endl;
		output << "///  Graft? " << std::boolalpha << options->design() << std::endl;
		output << "////  Weight: "<< options->weight() << std::endl;
		//output << std::endl;
		output << "///  SeqDesign? " << std::boolalpha << seq_options->design() << std::endl;
		if (seq_options->design()){
			output<< "////  Design Strategy: " << strategies[ seq_options->design_strategy() ] << std::endl;
		}
		if (! options->design()) continue;

		if (options->neighbor_min().size() >=1){
			utility::vector1<CDRNameEnum> neighbor_cdrs = options->neighbor_min();
			utility::vector1<std::string> neighbor_str;
			for (core::Size i = 1; i <= neighbor_cdrs.size(); ++i){
				neighbor_str.push_back(ab_info_->get_CDR_name(cdr));
			}
			output<< "///  Min Neighbors:  " << utility::to_string(neighbor_str) << std::endl;
		}

		output<< "///  Current Clusters only? " << std::boolalpha << set_options->include_only_current_cluster() << std::endl;
		output<< "///  Center Clusters only? " << std::boolalpha << set_options->include_only_center_clusters() << std::endl;
		if (cdr!=h3){
			output << "////" << std::endl;
			output << "////// Length Types: " << std::endl;
			output << "////" << std::endl;
			output << "///  1 " << std::boolalpha << set_options->length_type()[1] << std::endl;
			output << "///  2 " << std::boolalpha << set_options->length_type()[2] << std::endl;
			output << "///  3 " << std::boolalpha << set_options->length_type()[3] << std::endl;
		}
		output << "////" << std::endl;
		output << "////// Lengths: " << std::endl;
		output << "////" << std::endl;
		output << "/// Min " << set_options->min_length() << std::endl;
		output << "/// Max " << set_options->max_length() << std::endl;
		output << "///" << std::endl;

		//Include/Leave out options:

		vector1<std::string> exclude_clusters;
		vector1<std::string> include_clusters;
		vector1<CDRClusterEnum> exclude_clusters_enum = set_options->exclude_clusters();
		vector1<CDRClusterEnum> include_clusters_enum = set_options->include_only_clusters();

		for (core::Size i = 1; i <= exclude_clusters_enum.size(); ++i){
			exclude_clusters.push_back(ab_info_->get_cluster_name(exclude_clusters_enum[i]));
		}
		for (core::Size i = 1; i <= include_clusters_enum.size(); ++i){
			include_clusters.push_back(ab_info_->get_cluster_name(include_clusters_enum[i]));
		}

		output << "////// Include Only: " << std::endl;
		print_str_vec("PDBs", set_options->include_only_pdbs(), output);
		print_str_vec("Clusters", include_clusters);
		print_str_vec("Species", set_options->include_only_species(), output);
		print_str_vec("Germlines", set_options->include_only_germlines(), output);

		output<<"////// Exclude: " << std::endl;
		print_str_vec("PDBs", set_options->exclude_pdbs(), output);
		print_str_vec("Clusters", exclude_clusters, output);
		print_str_vec("Species", set_options->exclude_species(), output);
		print_str_vec("Germlines", set_options->exclude_germlines(), output);
		output<< std::endl;
	}
	output << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
}

void
AntibodyGraftDesignMover::print_str_vec(std::string const name, utility::vector1<std::string> const & vec, std::ostream & output) const {
	if (vec.size()>=1){
		output<< "/// "<< name <<" "<<utility::to_string(vec)<<std::endl;
	}
}

void
AntibodyGraftDesignMover::apply(pose::Pose & pose){

	if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
		utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	}
	//Create Instances.  Make sure everything is ready to begin.
	setup_scorefxn();
	setup_native_clusters(pose);
	setup_options_classes();
	seq_design_creator_ = AntibodySeqDesignTFCreatorOP( new AntibodySeqDesignTFCreator(ab_info_, cdr_seq_design_options_, 2 /* cdr_stem_size */) );

	if (cdr_set_.empty()){
		initialize_cdr_set(pose);
	}
	//Reinitialize.  Need to figure out how to do this properly via JD.
	top_designs_.clear();
	top_scores_.clear();

	//List of CDRs to graft.
	vector1< CDRNameEnum > cdrs_to_design;
	for (core::Size i=1; i<=CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		if (cdr_graft_design_options_[cdr]->design() ){
			if (cdr_set_[cdr].size() != 0){
				cdrs_to_design.push_back(cdr);
			}
			else {
				//Exit so that we don't continue thinking we are designing some cdr when we really arnt...
				utility_exit_with_message("CDR "  + ab_info_->get_CDR_name(cdr) + " set to graft design"
						" but no CDRs are in the set.  Please double check instructions settings.");
			}
		}
	}
	show(TR);
	if (cdrs_to_design.size() == 0){
		TR << "All CDRs fixed for low res graft designer or there are no CDRs in the set...." << std::endl;
		//Make sure our top designs has something in it for the antibody designer.
		check_for_top_designs(pose);
		return;
	}

	graft_mover_ = protocols::grafting::CCDEndsGraftMoverOP( new CCDEndsGraftMover(ab_info_->get_CDR_start(cdrs_to_design[1], pose)-1, ab_info_->get_CDR_end(cdrs_to_design[1], pose)+1) );
	anchored_graft_mover_ = protocols::grafting::AnchoredGraftMoverOP( new AnchoredGraftMover(ab_info_->get_CDR_start(cdrs_to_design[1], pose) - 1, ab_info_->get_CDR_end(cdrs_to_design[1], pose)+1) );
	setup_default_graft_settings();

	modeler_->set_scorefunction(scorefxn_);
	mover_generator_->set_scorefunction(scorefxn_);
	std::string ab_dock_chains = ab_info_->get_antibody_chain_string() +"_A";
	modeler_->ab_dock_chains(ab_dock_chains);
	mover_generator_->ab_dock_chains(ab_dock_chains);

	total_permutations_ = 1;
	for (core::Size i=1; i<=cdrs_to_design.size(); ++i){
		total_permutations_ *= cdr_set_[cdrs_to_design[i]].size();
	}
	TR<< "///// Total CDRs in set /////"<<std::endl;
	for(core::Size i = 1; i<=6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		modeler_->cdr_overhang(cdr, 2); // 2 Residue overhang for CDR packing and CDR minimization to go with grafting.  North CDR definitions already go into the beta barrel by a few residues
		mover_generator_->stem_size(2);
		TR << "/// "<<ab_info_->get_CDR_name(cdr)<<" "<<cdr_set_[cdr].size()<<std::endl;
	}

	TR <<"Total possible CDR combinations: "<< total_permutations_ << std::endl;

	paratope_epitope_cst_mover_ = constraints::ParatopeEpitopeSiteConstraintMoverOP( new ParatopeEpitopeSiteConstraintMover(ab_info_) );
	paratope_cst_mover_ = constraints::ParatopeSiteConstraintMoverOP( new ParatopeSiteConstraintMover(ab_info_) );
	setup_paratope_epitope_constraints(pose);
	mc_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo(pose, *scorefxn_, 1.0) );

	core::Real native_score = (*scorefxn_)(pose);
	scorefxn_->show(pose);
	//Temporary quick fix.

	if (benchmark_){
		this->setup_random_start_pose(pose, cdrs_to_design);
	}

	if (cdrs_to_design.size() == 1 && use_deterministic_algorithm_){
		run_deterministic_graft_algorithm(pose, cdrs_to_design);
	}
	else{
		//run_random_graft_algorithm(pose, cdrs_to_design);
		run_basic_mc_algorithm(pose, cdrs_to_design);
	}
	TR << "Native Pose: " << native_score << std::endl;
	TR << "Final Pose: " << (*scorefxn_)(pose) << std::endl;
	for (core::Size i = 2; i<= top_scores_.size(); ++i){
		TR << "Top Ensemble " << i << " : " << (*scorefxn_)(*top_designs_[i]) << std::endl;
	}
}

////////////////////////////////////////////// Boiler Plate ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
AntibodyGraftDesignMover::get_name() const{
	return "AntibodyGraftDesignMover";
}


} //design
} //antibody
} //protocols
