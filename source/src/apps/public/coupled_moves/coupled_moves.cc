// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <devel/init.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/CoupledMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/residue_selector/ClashBasedRepackShellSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyz.functions.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "apps.coupled_moves" );

OPT_1GRP_KEY(Integer, coupled_moves, ntrials)
OPT_1GRP_KEY(Real, coupled_moves, mc_kt)
OPT_1GRP_KEY(Real, coupled_moves, mm_bend_weight)
OPT_1GRP_KEY(Boolean, coupled_moves, trajectory)
OPT_1GRP_KEY(Boolean, coupled_moves, trajectory_gz)
OPT_1GRP_KEY(Integer, coupled_moves, trajectory_stride)
OPT_1GRP_KEY(String, coupled_moves, trajectory_file)
OPT_1GRP_KEY(String, coupled_moves, output_fasta)
OPT_1GRP_KEY(String, coupled_moves, output_stats)
OPT_1GRP_KEY(Boolean, coupled_moves, ligand_mode)
OPT_1GRP_KEY(Boolean, coupled_moves, initial_repack)
OPT_1GRP_KEY(Boolean, coupled_moves, save_sequences)
OPT_1GRP_KEY(Real, coupled_moves, ligand_prob)
OPT_1GRP_KEY(Boolean, coupled_moves, fix_backbone)
OPT_1GRP_KEY(Boolean, coupled_moves, bias_sampling)

void *
my_main( void* );

int
main( int argc, char * argv [] )
{
	try {

	OPT(in::path::database);
	OPT(in::ignore_unrecognized_res);
	OPT(out::nstruct);
	OPT(packing::resfile);
	OPT(in::file::native);
	NEW_OPT(coupled_moves::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(coupled_moves::mc_kt, "value of kT for Monte Carlo", 0.6);
	NEW_OPT(coupled_moves::mm_bend_weight, "weight of mm_bend bond angle energy term", 1.0);
	NEW_OPT(coupled_moves::trajectory, "record a trajectory", false);
	NEW_OPT(coupled_moves::trajectory_gz, "gzip the trajectory", false);
	NEW_OPT(coupled_moves::trajectory_stride, "write out a trajectory frame every N steps", 100);
	NEW_OPT(coupled_moves::trajectory_file, "name of trajectory file", "traj.pdb");
	NEW_OPT(coupled_moves::output_fasta, "name of FASTA output file", "sequences.fasta");
	NEW_OPT(coupled_moves::output_stats, "name of stats output file", "sequences.stats");
	NEW_OPT(coupled_moves::ligand_mode, "if true, model protein ligand interaction", false);
	NEW_OPT(coupled_moves::initial_repack, "start simulation with repack and design step", true);
	NEW_OPT(coupled_moves::save_sequences, "save all unique sequences", true);
	NEW_OPT(coupled_moves::ligand_prob, "probability of making a ligand move", 0.1);
	NEW_OPT(coupled_moves::fix_backbone, "do not make any backbone moves", false);
	NEW_OPT(coupled_moves::bias_sampling, "if true, bias rotamer selection based on energy", true);
	
	// initialize Rosetta
	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

	 } catch ( utility::excn::EXCN_Base const & e ) { 
		 std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

class CoupledMovesProtocol : public protocols::moves::Mover {
public:
	CoupledMovesProtocol();
	CoupledMovesProtocol(CoupledMovesProtocol const & cmp);

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "CoupledMovesProtocol"; }

	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new CoupledMovesProtocol( *this ) );
	}

	virtual protocols::moves::MoverOP	fresh_instance() const {
		return protocols::moves::MoverOP( new CoupledMovesProtocol );
	}

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	core::pack::task::TaskFactoryOP main_task_factory_;
};

CoupledMovesProtocol::CoupledMovesProtocol(): Mover(),
		score_fxn_(core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() )),
		main_task_factory_(core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory() ))
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	main_task_factory_->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( TaskOperationCOP( new operation::ReadResfile ) );
	} else {
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory_->push_back( rtrop );
	}
	
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	//main_task_factory_->push_back( new operation::PreserveCBeta );

	// set up the score function and add the bond angle energy term
	score_fxn_ = core::scoring::get_score_function();
	score_fxn_->set_weight(core::scoring::mm_bend, option[ coupled_moves::mm_bend_weight ]);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn_->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ backrub::pivot_atoms ]);
	score_fxn_->set_energy_method_options(energymethodoptions);
}

CoupledMovesProtocol::CoupledMovesProtocol(CoupledMovesProtocol const & cmp): Mover(cmp),
		score_fxn_(cmp.score_fxn_),
		main_task_factory_(cmp.main_task_factory_)
{}

void CoupledMovesProtocol::apply( core::pose::Pose& pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	score_fxn_->show(TR, pose);
	TR.flush();

	protocols::moves::MonteCarlo mc(pose, *score_fxn_, option[ coupled_moves::mc_kt ]);

	protocols::viewer::add_monte_carlo_viewer(mc, "CoupledMoves", 600, 600);
	
	////////////////// Material above is scheduled for constructor //////////////////

	std::string output_tag(protocols::jd2::current_output_name());

	// start with a fresh copy of the optimized pose
	core::pose::PoseOP pose_copy(new core::pose::Pose(pose));

	// reset the Monte Carlo object
	mc.reset(*pose_copy);
	
	core::pack::task::PackerTaskOP task( main_task_factory_->create_task_and_apply_taskoperations( *pose_copy ) );
	
	protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
	if (option[ coupled_moves::trajectory ]) {
		trajectory.file_name(output_tag + "_traj.pdb" + (option[ coupled_moves::trajectory_gz ] ? ".gz" : ""));
		trajectory.stride(option[ coupled_moves::trajectory_stride ]);
		trajectory.reset(mc);
	}

	utility::vector1<core::Size> move_positions;
	utility::vector1<core::Size> design_positions;
	
	// using ClashBasedRepackShellSelector to define repack shell
	core::pack::task::residue_selector::ClashBasedRepackShellSelectorOP rs( new core::pack::task::residue_selector::ClashBasedRepackShellSelector(task, score_fxn_));
	utility::vector1< bool > to_repack = rs->apply( *pose_copy );

	for(core::Size i = 1; i <= to_repack.size(); ++i) {
		if (task->residue_task(i).has_behavior("AUTO")) {
			if (to_repack[i]) {
				task->nonconst_residue_task(i).restrict_to_repacking();
			} else {
				task->nonconst_residue_task(i).prevent_repacking();
			}
		}
	}
	
	for(core::Size i = 1; i <= to_repack.size(); ++i) {
		if (task->design_residue(i)) {
			design_positions.push_back(i);
			move_positions.push_back(i);
		}
		else if (task->pack_residue(i)) {
			move_positions.push_back(i);
		}
	}
	
	TR << "Design positions ";
	for(core::Size i = 1; i <= design_positions.size(); i++) {
		TR << design_positions[i] << " ";
	}
	TR << std::endl;
	
	// ASSUMPTION: the ligand is the last residue in the given PDB
	core::Size lig_res = pose_copy->total_residue();
	
	protocols::simple_moves::CoupledMoverOP coupled_mover;
	
	if ( option[coupled_moves::ligand_mode] ) {
		coupled_mover = protocols::simple_moves::CoupledMoverOP(new protocols::simple_moves::CoupledMover(pose_copy, score_fxn_, task, lig_res));
	} else {
		coupled_mover = protocols::simple_moves::CoupledMoverOP(new protocols::simple_moves::CoupledMover(pose_copy, score_fxn_, task));
	}
	
	coupled_mover->set_fix_backbone( option[coupled_moves::fix_backbone] );
	coupled_mover->set_bias_sampling( option[coupled_moves::bias_sampling] );
	
	protocols::simple_moves::PackRotamersMoverOP pack(new protocols::simple_moves::PackRotamersMover( score_fxn_, task, 1 ));		

	if ( option[coupled_moves::initial_repack] ) {
		pack->apply(*pose_copy);
	}
	
	TR << "Running " << option[ coupled_moves::ntrials ] << " trials..." << std::endl;

	std::map<std::string,core::Real> unique_sequences;
	std::map<std::string,core::scoring::EnergyMap> unique_scores;
	
	core::Size ntrials = option[ coupled_moves::ntrials ];
	
	std::string resfile_name = option[packing::resfile]()[1];
	
	for (core::Size i = 1; i <= ntrials; ++i) {
		core::Size random = numeric::random::random_range(1, move_positions.size());
		core::Size resnum = move_positions[random];
		std::string move_type;
		core::Real move_prob = numeric::random::uniform();
		if (move_prob < option[coupled_moves::ligand_prob]) {
			resnum = lig_res;
			move_type = "LIGAND";
		} else {
			move_type = "RESIDUE";
		}
		
		coupled_mover->set_resnum(resnum);
		coupled_mover->apply(*pose_copy);
		
		bool accepted = mc.boltzmann(*pose_copy, move_type);
		
		if (accepted) {
			core::Real current_score = mc.last_accepted_score();
			std::string sequence = "";
			for(core::Size index = 1; index <= design_positions.size(); index++) {
				sequence += pose_copy->residue(design_positions[index]).name1();
			}
			
			TR << i << " " << sequence << " " << pose_copy->energies().total_energy() << std::endl;
			
			if (option[coupled_moves::save_sequences]) {
				if (unique_sequences.find(sequence) == unique_sequences.end()) {
					unique_sequences.insert(std::make_pair(sequence,current_score));
					unique_scores.insert(std::make_pair(sequence,pose_copy->energies().total_energies()));
				} else {
					if (unique_sequences[sequence] > current_score) {
						unique_sequences[sequence] = current_score;
						unique_scores[sequence] = pose_copy->energies().total_energies();
					}
				}
			}
			
		}
				
		if (option[ coupled_moves::trajectory ]) trajectory.update_after_boltzmann(mc);
		
	}
	
	mc.show_counters();

	// dump out the low score and last accepted poses
	
	TR << "Last Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	*pose_copy = mc.lowest_score_pose();

	TR << "Low Score:" << std::endl;
	score_fxn_->show(TR, *pose_copy);
	TR.flush();

	pose = mc.lowest_score_pose();

	mc.last_accepted_pose().dump_pdb(output_tag + "_last.pdb");
	mc.lowest_score_pose().dump_pdb(output_tag + "_low.pdb");

	if (option[coupled_moves::save_sequences]) {
		std::ofstream out_fasta( (output_tag + ".fasta").c_str() );
		core::Size count = 1;
		for(std::map<std::string,core::Real>::iterator it = unique_sequences.begin(); it != unique_sequences.end(); it++) {
			out_fasta << ">Sequence" << count << " " << it->second << std::endl;
			out_fasta << it->first << std::endl;
			count++;
		}
		out_fasta.close();
			
		std::ofstream out_stats( (output_tag + ".stats").c_str() );
		count = 1;
		for(std::map<std::string,core::Real>::iterator it = unique_sequences.begin(); it != unique_sequences.end(); it++) {
			out_stats << "Sequence" << count << "\t" << it->second << "\tsequence:\t" << it->first << "\t" << unique_scores[it->first].weighted_string_of(score_fxn_->weights()) << std::endl;
			count++;
		}
		out_stats.close();
	}

}

typedef utility::pointer::shared_ptr<CoupledMovesProtocol> CoupledMovesProtocolOP;

void *
my_main( void* )
{
	
	CoupledMovesProtocolOP coupled_moves( new CoupledMovesProtocol );
	protocols::jd2::JobDistributor::get_instance()->go( coupled_moves );
	
	return 0;
}
