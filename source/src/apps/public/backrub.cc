// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file backrub.cc
/// @brief run backrub Monte Carlo
/// @author Colin A. Smith (colin.smith@ucsf.edu)
/// @detailed
/// Currently a work in progress. The goal is to match the features of rosetta++ -backrub_mc

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/keys/Key3Vector.hh>

static numeric::random::RandomGenerator RG(222578262);
basic::Tracer TR("apps.backrub");

OPT_1GRP_KEY(Integer, backrub, ntrials)
OPT_1GRP_KEY(Real, backrub, sc_prob)
OPT_1GRP_KEY(Real, backrub, sm_prob)
OPT_1GRP_KEY(Real, backrub, sc_prob_uniform)
OPT_1GRP_KEY(Real, backrub, sc_prob_withinrot)
OPT_1GRP_KEY(Real, backrub, mc_kt)
OPT_1GRP_KEY(Real, backrub, mm_bend_weight)
OPT_1GRP_KEY(Boolean, backrub, initial_pack)
OPT_1GRP_KEY(File, backrub, minimize_movemap)
OPT_1GRP_KEY(Boolean, backrub, trajectory)
OPT_1GRP_KEY(Boolean, backrub, trajectory_gz)
OPT_1GRP_KEY(Integer, backrub, trajectory_stride)

void *
my_main( void* );

int
main( int argc, char * argv [] )
{
	try {

	OPT(in::path::database);
	OPT(in::file::s);
	OPT(in::file::l);
	OPT(in::file::movemap);
	OPT(in::ignore_unrecognized_res);
	OPT(out::nstruct);
	OPT(packing::resfile);
	OPT(constraints::cst_fa_weight);
	OPT(constraints::cst_fa_file);
	OPT(backrub::pivot_residues);
	OPT(backrub::pivot_atoms);
	OPT(backrub::min_atoms);
	OPT(backrub::max_atoms);
	NEW_OPT(backrub::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(backrub::sc_prob, "probability of making a side chain move", 0.25);
	NEW_OPT(backrub::sm_prob, "probability of making a small move", 0);
	NEW_OPT(backrub::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.1);
	NEW_OPT(backrub::sc_prob_withinrot, "probability of sampling within the current rotamer", 0.0);
	NEW_OPT(backrub::mc_kt, "value of kT for Monte Carlo", 0.6);
	NEW_OPT(backrub::mm_bend_weight, "weight of mm_bend bond angle energy term", 1.0);
	NEW_OPT(backrub::initial_pack, "force a repack at the beginning regardless of whether mutations are set in the resfile", false);
	NEW_OPT(backrub::minimize_movemap, "specify degrees of freedom for minimization", "");
	NEW_OPT(backrub::trajectory, "record a trajectory", false);
	NEW_OPT(backrub::trajectory_gz, "gzip the trajectory", false);
	NEW_OPT(backrub::trajectory_stride, "write out a trajectory frame every N steps", 100);

	// initialize Rosetta
	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

	 } catch ( utility::excn::EXCN_Base const & e ) { 
		 std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

bool
read_fold_tree_from_file(
	core::kinematics::FoldTree & foldtree,
	std::string filepath
)
{
	std::ifstream filestream(filepath.c_str());

	while (filestream.good()) {

		std::string line;
		std::string key;

		getline(filestream, line);
		if (filestream.fail()) {
			//TR << "getline() failed" << std::endl;
			return false;
		}

		std::istringstream linestream(line);
		linestream >> key;
		if (key == "FOLD_TREE") {
			linestream.clear();
			linestream.seekg(0, std::ios::beg);
			linestream >> foldtree;
			if (linestream.fail()) {
				TR << "FoldTree parsing failed" << std::endl;
				return false;
			} else {
				return true;
			}
		}
	}

	return false;
}

bool
read_fold_tree_from_file(
	core::pose::Pose & pose,
	std::string filepath
)
{
	core::kinematics::FoldTree foldtree;

	if (read_fold_tree_from_file(foldtree, filepath)) {
		if (foldtree.nres() == pose.total_residue()) {
			pose.fold_tree(foldtree);
			return true;
		} else {
			TR << "Different number of residues in Pose (" << pose.total_residue() << ") and FoldTree (" << foldtree.nres()
			   << ")" << std::endl;
		}
	}

	return false;
}

void
append_fold_tree_to_file(
	core::kinematics::FoldTree const & foldtree,
	std::string file_path
)
{
	std::ofstream filestream(file_path.c_str(), std::ios::out|std::ios::app);
	if (filestream.good()) {
		filestream << foldtree;
		filestream.close();
	} else {
		TR << "couldn't open file to append FoldTree" << std::endl;
	}
}

utility::vector1<core::Size>
positions_incompatible_with_task(
	core::pose::Pose & pose,
	core::pack::task::PackerTask & packertask
)
{
	utility::vector1<core::Size> incompatible_positions;

	assert(pose.total_residue() == packertask.total_residue());

	// iterate over all residues to see if they're compatible
	for (core::Size i = 1; i <= pose.total_residue(); ++i) {

		// only check packable residues for compatibility
		if (packertask.pack_residue(i)) {

			// assume residue is incompatible
			bool incompatible(true);

			// check to see if pose residue type is in list of allowed residue types
			core::pack::task::ResidueLevelTask const & residueleveltask(packertask.residue_task(i));
			for (core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter iter(residueleveltask.allowed_residue_types_begin());
			     iter != residueleveltask.allowed_residue_types_end(); ++iter) {

				if ((*iter)->name() == pose.residue_type(i).name()) incompatible = false;
			}

			if (incompatible) incompatible_positions.push_back(i);
		}
	}

	return incompatible_positions;
}

class BackrubProtocol : public protocols::moves::Mover {
public:
	BackrubProtocol();
	BackrubProtocol(BackrubProtocol const & bp);

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "BackrubProtocol"; }

	virtual protocols::moves::MoverOP clone() const {
		return new BackrubProtocol( *this );
	}

	virtual	 protocols::moves::MoverOP	fresh_instance() const {
		return new BackrubProtocol;
	}

	protocols::backrub::BackrubMoverOP get_backrub_mover() const;

	void write_database(){
		backrubmover_->branchopt().write_database();
	}

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;
	protocols::backrub::BackrubMoverOP backrubmover_;
	protocols::simple_moves::SmallMover smallmover_;
	protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover_;
	protocols::simple_moves::PackRotamersMover packrotamersmover_;
};

BackrubProtocol::BackrubProtocol(): Mover(),
		score_fxn_(new core::scoring::ScoreFunction()),
		main_task_factory_(new core::pack::task::TaskFactory()),
		backrubmover_(new protocols::backrub::BackrubMover()),
		smallmover_(protocols::simple_moves::SmallMover()),
		sidechainmover_(protocols::simple_moves::sidechain_moves::SidechainMover()),
		packrotamersmover_(protocols::simple_moves::PackRotamersMover())
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pack::task;

	main_task_factory_->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( new operation::ReadResfile );
	} else {
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory_->push_back( rtrop );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory_->push_back( new operation::PreserveCBeta );

	// set up the score function and add the bond angle energy term
	score_fxn_ = core::scoring::getScoreFunction();
	score_fxn_->set_weight(core::scoring::mm_bend, option[ backrub::mm_bend_weight ]);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn_->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ backrub::pivot_atoms ]);
	score_fxn_->set_energy_method_options(energymethodoptions);
	if ( option[ in::file::centroid_input ].user() ) {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	} else {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	}

	// read known and unknown optimization parameters from the database
	backrubmover_->branchopt().read_database();
	// tell the branch angle optimizer about the score function MMBondAngleResidueTypeParamSet, if any
	if (energymethodoptions.bond_angle_residue_type_param_set()) {
		backrubmover_->branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
	}

	// set up the SmallMover
	smallmover_.nmoves(1);
	if (option[ backrub::sm_prob ] > 0) {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->init_from_file(option[ in::file::movemap ]);
		smallmover_.movemap(movemap);
	}

	// set up the SidechainMover
	sidechainmover_.set_task_factory(main_task_factory_);
	sidechainmover_.set_prob_uniform(option[ backrub::sc_prob_uniform ]);
	sidechainmover_.set_prob_withinrot(option[ backrub::sc_prob_withinrot ]);

	// set up the PackRotamersMover
	packrotamersmover_.task_factory(main_task_factory_);
	packrotamersmover_.score_function(score_fxn_);

	//backrubmover_= new protocols::backrub::BackrubMover;
}

BackrubProtocol::BackrubProtocol(BackrubProtocol const & bp): Mover(bp),
		score_fxn_(bp.score_fxn_),
		main_task_factory_(bp.main_task_factory_),
		backrubmover_(bp.backrubmover_),
		smallmover_(bp.smallmover_),
		sidechainmover_(bp.sidechainmover_),
		packrotamersmover_(bp.packrotamersmover_)
{}


protocols::backrub::BackrubMoverOP BackrubProtocol::get_backrub_mover() const{
	return backrubmover_;
}

void BackrubProtocol::apply( core::pose::Pose& pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool custom_fold_tree = false;
	std::string input_tag = protocols::jd2::JobDistributor::get_instance()->current_job()->inner_job()->input_tag();
	if ( option[ in::file::centroid_input ].user() ) {
		TR.Warning << "*** This is untested with centroid mode! ***" << std::endl;
		//core::import_pose::centroid_pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() );
		custom_fold_tree = read_fold_tree_from_file( pose, input_tag );
		core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
	} else {
		//core::import_pose::pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() );
		custom_fold_tree = read_fold_tree_from_file( pose, input_tag );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);
	}

	backrubmover_->clear_segments();
	backrubmover_->set_input_pose( pose.clone() ); // hmm, is that okay? NO NO NO it is not.
	backrubmover_->add_mainchain_segments_from_options();

	TR << "Score After PDB Load:" << std::endl;
	score_fxn_->show(TR, pose);
	TR.flush();

	backrubmover_->optimize_branch_angles(pose);
	//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postoptbranch.pdb");
	sidechainmover_.idealize_sidechains(pose);
	//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postidealizesc.pdb");
	core::Real sc_prob = option[ backrub::sc_prob ];
	if (sc_prob && sidechainmover_.packed_residues().size() == 0) {
		sc_prob = 0;
		TR << "Warning: No side chains to move. Not using SidechainMover." << std::endl;
	}

	TR << "Score After Branch Angle Optimization/Side Chain Idealization:" << std::endl;
	score_fxn_->show(TR, pose);
	TR.flush();

	// check to see if we need to force a repack
	bool initial_pack(option[ backrub::initial_pack ]);
	core::pack::task::PackerTaskOP temp_task(main_task_factory_->create_task_and_apply_taskoperations(pose));
	utility::vector1<core::Size> incompatible_positions(positions_incompatible_with_task(pose, *temp_task));
	if (incompatible_positions.size()) {
		TR << "Starting ResidueType not allowed in resfile at position(s):";
		for (core::Size i = 1; i <= incompatible_positions.size(); i++) {
			TR << " " << incompatible_positions[i];
		}
		TR << std::endl;
		if (!initial_pack) {
			initial_pack = true;
			TR << "Forcing initial pack" << std::endl;
		}
	}

	protocols::moves::MonteCarlo mc(pose, *score_fxn_, option[ backrub::mc_kt ]);

	// create viewer windows if OpenGL is enabled
	protocols::viewer::add_monte_carlo_viewer(mc, "Backrub", 600, 600);

	////////////////// Material above is scheduled for constructor //////////////////


	std::string output_tag(protocols::jd2::current_output_name());

	// start with a fresh copy of the optimized pose
	core::pose::PoseOP pose_copy(new core::pose::Pose(pose));

	// repack/redesign at the beginning if specified/necessary
	if (initial_pack) {
		packrotamersmover_.apply(*pose_copy);
		//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postpack.pdb");

		// if a minimization movemap was specified, go through a series of minimizations
		if (option[ backrub::minimize_movemap ].user()) {

			// setup the MoveMaps
			core::kinematics::MoveMapOP minimize_movemap = new core::kinematics::MoveMap;
			minimize_movemap->init_from_file(option[ backrub::minimize_movemap ]);
			core::kinematics::MoveMapOP minimize_movemap_progressive = new core::kinematics::MoveMap;

			// setup the MinMover
			protocols::simple_moves::MinMover minmover;
			minmover.score_function(score_fxn_);
			minmover.min_type("dfpmin");

			// first minimize just the side chains
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
			     iter != minimize_movemap->movemap_torsion_id_end(); ++iter) {
				if (iter->first.second == core::id::CHI) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminchi.pdb");

			// next minimize the side chains and backbone
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
			     iter != minimize_movemap->movemap_torsion_id_end(); ++iter) {
				if (iter->first.second == core::id::BB) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminbb.pdb");

			// finally minimize everything
			for (core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
			     iter != minimize_movemap->movemap_torsion_id_end(); ++iter) {
				if (iter->first.second == core::id::JUMP) minimize_movemap_progressive->set(iter->first, iter->second);
			}
			minmover.movemap(minimize_movemap_progressive);
			minmover.apply(*pose_copy);
			//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminjump.pdb");
		}
	}



	// reset the Monte Carlo object
	mc.reset(*pose_copy);

	protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
	if (option[ backrub::trajectory ]) {
		trajectory.file_name(output_tag + "_traj.pdb" + (option[ backrub::trajectory_gz ] ? ".gz" : ""));
		trajectory.stride(option[ backrub::trajectory_stride ]);
		trajectory.reset(mc);
	}

	TR << "Running " << option[ backrub::ntrials ] << " trials..." << std::endl;

	for (int i = 1; i <= option[ backrub::ntrials ]; ++i) {

		std::string move_type;

		// could use random mover for this...
		core::Real move_prob = RG.uniform();
		if (move_prob > option[ backrub::sm_prob ] + option[ backrub::sc_prob ]) {
			backrubmover_->apply(*pose_copy);
			move_type = backrubmover_->type();
		} else if (move_prob > sc_prob) {
			smallmover_.apply(*pose_copy);
			move_type = smallmover_.type();
		} else {
			sidechainmover_.apply(*pose_copy);
			move_type = sidechainmover_.type();
		}

		mc.boltzmann(*pose_copy, move_type);

		if (option[ backrub::trajectory ]) trajectory.update_after_boltzmann(mc);
	}

	mc.show_counters();

	// repack II: LAST
	// this could also be done IN the loop, so that we get a minimized structure right after the monte carlo run
	//packrotamersmover.apply(*pose);

	// dump out the low score and last accepted poses
	TR << "Last Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	*pose_copy = mc.lowest_score_pose();

	// repack II: LOW
	//packrotamersmover.apply(*pose);

	TR << "Low Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	pose= mc.lowest_score_pose();

	mc.last_accepted_pose().dump_pdb(output_tag + "_last.pdb");
	mc.lowest_score_pose().dump_pdb(output_tag + "_low.pdb");

	if (custom_fold_tree) {
		append_fold_tree_to_file(mc.lowest_score_pose().fold_tree(),  output_tag + "_low.pdb"); // this is the lowest scoring from MC trials
		append_fold_tree_to_file(mc.last_accepted_pose().fold_tree(), output_tag + "_last.pdb"); // this is the last accepted pose from MC trials
	}
}

typedef utility::pointer::owning_ptr<BackrubProtocol> BackrubProtocolOP;

void *
my_main( void* )
{

	BackrubProtocolOP backrub_protocol = new BackrubProtocol;
	protocols::jd2::JobDistributor::get_instance()->go( backrub_protocol );

	// write parameters for any sets of branching atoms for which there were not optimization coefficients
	backrub_protocol->write_database();
	return 0;
}
