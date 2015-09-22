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
/// @details
/// Currently a work in progress. The goal is to match the features of rosetta++ -backrub_mc

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/DOFHistogramRecorder.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/viewer/viewers.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/DOF_ID_Range.hh>
#include <devel/init.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
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
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

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


static THREAD_LOCAL basic::Tracer TR( "apps.backrub" );

//JAB - This is a fork of Backrub, and looks like it acted as a test for the detailed balance paper and code.
// It should probably be deprecated.

OPT_1GRP_KEY(Real, backrub, backrub_sc_prob)
OPT_1GRP_KEY(Real, backrub, sm_angle_max)
OPT_1GRP_KEY(Boolean, backrub, detailed_balance)
OPT_1GRP_KEY(Boolean, backrub, test_detailed_balance)


void *
my_main( void* );

int
main( int argc, char * argv [] )
{ try {
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
	NEW_OPT(backrub::backrub_sc_prob, "probability of making a side chain move during 3 residue backrub moves", 0);
	NEW_OPT(backrub::sm_angle_max, "small move maximum band of angluar perturbation", 6);
	NEW_OPT(backrub::detailed_balance, "preserve detailed balance", false);
	NEW_OPT(backrub::test_detailed_balance, "test that detailed balance is preserved", false);

	// initialize Rosetta
	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
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

	while ( filestream.good() ) {

		std::string line;
		std::string key;

		getline(filestream, line);
		if ( filestream.fail() ) {
			//TR << "getline() failed" << std::endl;
			return false;
		}

		std::istringstream linestream(line);
		linestream >> key;
		if ( key == "FOLD_TREE" ) {
			linestream.clear();
			linestream.seekg(0, std::ios::beg);
			linestream >> foldtree;
			if ( linestream.fail() ) {
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

	if ( read_fold_tree_from_file(foldtree, filepath) ) {
		if ( foldtree.nres() == pose.total_residue() ) {
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
	if ( filestream.good() ) {
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
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

		// only check packable residues for compatibility
		if ( packertask.pack_residue(i) ) {

			// assume residue is incompatible
			bool incompatible(true);

			// check to see if pose residue type is in list of allowed residue types
			core::pack::task::ResidueLevelTask const & residueleveltask(packertask.residue_task(i));
			for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter iter(residueleveltask.allowed_residue_types_begin());
					iter != residueleveltask.allowed_residue_types_end(); ++iter ) {

				if ( (*iter)->name() == pose.residue_type(i).name() ) incompatible = false;
			}

			if ( incompatible ) incompatible_positions.push_back(i);
		}
	}

	return incompatible_positions;
}

void *
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// create a TaskFactory with the resfile
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	TaskFactoryOP main_task_factory( new TaskFactory );
	main_task_factory->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( TaskOperationCOP( new operation::ReadResfile ) );
	} else {
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory->push_back( rtrop );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( TaskOperationCOP( new operation::PreserveCBeta ) );

	// set up the score function and add the bond angle energy term
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
	score_fxn->set_weight(core::scoring::mm_bend, option[ backrub::mm_bend_weight ]);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn->energy_method_options());
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ backrub::pivot_atoms ]);
	score_fxn->set_energy_method_options(energymethodoptions);
	if ( option[ in::file::centroid_input ].user() ) {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn);
	} else {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	}
	// use an empty score function if testing detailed balance
	if ( option[ backrub::test_detailed_balance ] ) score_fxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );

	// set up the BackrubMover
	protocols::backrub::BackrubMover backrubmover;
	// read known and unknown optimization parameters from the database
	backrubmover.branchopt().read_database();
	// tell the branch angle optimizer about the score function MMBondAngleResidueTypeParamSet, if any
	if ( energymethodoptions.bond_angle_residue_type_param_set() ) {
		backrubmover.branchopt().bond_angle_residue_type_param_set(energymethodoptions.bond_angle_residue_type_param_set());
	}
	backrubmover.set_preserve_detailed_balance(option[ backrub::detailed_balance ]);

	// set up the SmallMover
	protocols::simple_moves::SmallMover smallmover;
	smallmover.nmoves(1);
	smallmover.set_preserve_detailed_balance(option[ backrub::detailed_balance ]);
	if ( option[ backrub::sm_angle_max ].user() ) smallmover.angle_max(option[ backrub::sm_angle_max ]);
	if ( option[ backrub::sm_prob ] > 0 ) {
		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
		movemap->init_from_file(option[ in::file::movemap ]);
		smallmover.movemap(movemap);
	}

	// set up the SidechainMover
	protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover;
	sidechainmover.set_task_factory(main_task_factory);
	sidechainmover.set_prob_uniform(option[ backrub::sc_prob_uniform ]);
	sidechainmover.set_prob_withinrot(option[ backrub::sc_prob_withinrot ]);
	sidechainmover.set_preserve_detailed_balance(option[ backrub::detailed_balance ]);

	// set up the PackRotamersMover
	protocols::simple_moves::PackRotamersMover packrotamersmover;
	packrotamersmover.task_factory(main_task_factory);
	packrotamersmover.score_function(score_fxn);

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	for ( core::Size jobnum = 1; jobnum <= input_jobs.size(); ++jobnum ) {

		TR << "Processing " << input_jobs[jobnum]->input_tag() << "..." << std::endl;
		bool custom_fold_tree(false);

		// load the PDB file
		core::pose::PoseOP input_pose( new core::pose::Pose() );
		if ( option[ in::file::centroid_input ].user() ) {
			TR.Warning << "*** This is untested with centroid mode! ***" << std::endl;
			core::import_pose::centroid_pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() );
			custom_fold_tree = read_fold_tree_from_file( *input_pose, input_jobs[jobnum]->input_tag() );
			core::scoring::constraints::add_constraints_from_cmdline_to_pose(*input_pose);
		} else {
			core::import_pose::pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() );
			custom_fold_tree = read_fold_tree_from_file( *input_pose, input_jobs[jobnum]->input_tag() );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose(*input_pose);
		}
		//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postread.pdb");

		backrubmover.clear_segments();
		backrubmover.set_input_pose(input_pose);
		backrubmover.add_mainchain_segments_from_options();

		TR << "Score After PDB Load:" << std::endl;
		score_fxn->show(TR, *input_pose);
		TR.flush();

		backrubmover.optimize_branch_angles(*input_pose);
		//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postoptbranch.pdb");
		sidechainmover.idealize_sidechains(*input_pose);
		//input_pose->dump_pdb(input_jobs[jobnum]->output_tag(0) + "_postidealizesc.pdb");
		core::Real sc_prob = option[ backrub::sc_prob ];
		if ( sc_prob && sidechainmover.packed_residues().size() == 0 ) {
			sc_prob = 0;
			TR << "Warning: No side chains to move. Not using SidechainMover." << std::endl;
		}

		TR << "Score After Branch Angle Optimization/Side Chain Idealization:" << std::endl;
		score_fxn->show(TR, *input_pose);
		TR.flush();

		// check to see if we need to force a repack
		bool initial_pack(option[ backrub::initial_pack ]);
		core::pack::task::PackerTaskOP temp_task(main_task_factory->create_task_and_apply_taskoperations(*input_pose));
		utility::vector1<core::Size> incompatible_positions(positions_incompatible_with_task(*input_pose, *temp_task));
		if ( incompatible_positions.size() ) {
			TR << "Starting ResidueType not allowed in resfile at position(s):";
			for ( core::Size i = 1; i <= incompatible_positions.size(); i++ ) {
				TR << " " << incompatible_positions[i];
			}
			TR << std::endl;
			if ( !initial_pack ) {
				initial_pack = true;
				TR << "Forcing initial pack" << std::endl;
			}
		}

		protocols::moves::MonteCarlo mc(*input_pose, *score_fxn, option[ backrub::mc_kt ]);

		// create viewer windows if OpenGL is enabled
		protocols::viewer::add_monte_carlo_viewer(mc, "Backrub", 600, 600);

		// iterate to generate multiple structures
		for ( int structnum = 1; structnum <= input_jobs[jobnum]->nstruct(); ++structnum ) {
			// start with a fresh copy of the optimized pose
			core::pose::PoseOP pose( new core::pose::Pose(*input_pose) );

			// repack/redesign at the beginning if specified/necessary
			if ( initial_pack ) {
				packrotamersmover.apply(*pose);
				//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postpack.pdb");

				// if a minimization movemap was specified, go through a series of minimizations
				if ( option[ backrub::minimize_movemap ].user() ) {

					// setup the MoveMaps
					core::kinematics::MoveMapOP minimize_movemap( new core::kinematics::MoveMap );
					minimize_movemap->init_from_file(option[ backrub::minimize_movemap ]);
					core::kinematics::MoveMapOP minimize_movemap_progressive( new core::kinematics::MoveMap );

					// setup the MinMover
					protocols::simple_moves::MinMover minmover;
					minmover.score_function(score_fxn);
					minmover.min_type("dfpmin");

					// first minimize just the side chains
					for ( core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
							iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
						if ( iter->first.second == core::id::CHI ) minimize_movemap_progressive->set(iter->first, iter->second);
					}
					minmover.movemap(minimize_movemap_progressive);
					minmover.apply(*pose);
					//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminchi.pdb");

					// next minimize the side chains and backbone
					for ( core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
							iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
						if ( iter->first.second == core::id::BB ) minimize_movemap_progressive->set(iter->first, iter->second);
					}
					minmover.movemap(minimize_movemap_progressive);
					minmover.apply(*pose);
					//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminbb.pdb");

					// finally minimize everything
					for ( core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator iter = minimize_movemap->movemap_torsion_id_begin();
							iter != minimize_movemap->movemap_torsion_id_end(); ++iter ) {
						if ( iter->first.second == core::id::JUMP ) minimize_movemap_progressive->set(iter->first, iter->second);
					}
					minmover.movemap(minimize_movemap_progressive);
					minmover.apply(*pose);
					//pose->dump_pdb(input_jobs[jobnum]->output_tag(structnum) + "_postminjump.pdb");
				}
			}

			std::string output_tag(input_jobs[jobnum]->output_tag(structnum));

			// reset the Monte Carlo object
			mc.reset(*pose);

			protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
			if ( option[ backrub::trajectory ] ) {
				trajectory.file_name(output_tag + "_traj.pdb" + (option[ backrub::trajectory_gz ] ? ".gz" : ""));
				trajectory.stride(option[ backrub::trajectory_stride ]);
				trajectory.reset(mc);
			}

			//numeric::MultiDimensionalHistogram mdhist(pose->residue(1).chi().size());
			//mdhist.set_dimension(1, 12, -180, 180, "chi1");
			//mdhist.set_dimension(2, 12, -180, 180, "chi2");
			//mdhist.set_dimension(3, 12, -180, 180, "chi3");
			//mdhist.set_dimension(4, 12, -180, 180, "chi4");
			//numeric::MultiDimensionalHistogram mdhist_proposed(mdhist);
			//mdhist.label("trajectory");
			//mdhist_proposed.label("proposed");

			protocols::simple_moves::DOFHistogramRecorder dof_recorder;

			if ( option[ backrub::test_detailed_balance ] ) {

				std::set<core::id::DOF_ID_Range> range_set;

				if ( option[ backrub::sm_prob ] + option[ backrub::sc_prob ] < 1 ) {
					utility::vector1<core::id::DOF_ID_Range> backrub_ranges(backrubmover.dof_id_ranges(*pose));
					range_set.insert(backrub_ranges.begin(), backrub_ranges.end());
				}

				if ( option[ backrub::sc_prob ] ) {
					utility::vector1<core::id::DOF_ID_Range> sidechainmover_ranges(sidechainmover.dof_id_ranges(*pose));
					range_set.insert(sidechainmover_ranges.begin(), sidechainmover_ranges.end());
				}

				if ( option[ backrub::sm_prob ] ) {
					utility::vector1<core::id::DOF_ID_Range> smallmover_ranges(smallmover.dof_id_ranges(*pose));
					range_set.insert(smallmover_ranges.begin(), smallmover_ranges.end());
				}

				utility::vector1<core::id::DOF_ID_Range> range_vector(range_set.begin(), range_set.end());

				dof_recorder.num_bins(12);
				dof_recorder.insert_dofs_by_residue(*pose, range_vector);

				//TR << "DOFs to Be Modified" << std::endl;
				//for (core::Size i = 1; i <= range_vector.size(); ++i) {
				// TR << range_vector[i] << std::endl;
				//}
			}

			TR << "Running " << option[ backrub::ntrials ] << " trials..." << std::endl;

			for ( int i = 1; i <= option[ backrub::ntrials ]; ++i ) {

				std::string move_type;
				core::Real proposal_density_ratio(1);

				// could use random mover for this...
				core::Real move_prob = numeric::random::rg().uniform();
				if ( move_prob > option[ backrub::sm_prob ] + option[ backrub::sc_prob ] ) {
					backrubmover.apply(*pose);
					move_type = backrubmover.type();

					protocols::backrub::BackrubSegment const & segment(backrubmover.segment(backrubmover.last_segment_id()));
					if ( option[ backrub::backrub_sc_prob ] && segment.size() == 7 ) {
						core::Size middle_resnum(static_cast<core::Size>((segment.start_atomid().rsd() + segment.end_atomid().rsd())*.5));
						//TR << "Simultaneous move for segment: " << segment.start_atomid() << segment.end_atomid() << " middle: " << middle_resnum << std::endl;
						if ( sidechainmover.residue_packed()[middle_resnum] ) {
							if ( numeric::random::rg().uniform() < option[ backrub::backrub_sc_prob ] ) {
								sidechainmover.next_resnum(middle_resnum);
								//TR << "next_resnum before: " << sidechainmover.next_resnum() << std::endl;
								sidechainmover.apply(*pose);
								//TR << "next_resnum after: " << sidechainmover.next_resnum() << std::endl;
								move_type += "_" + sidechainmover.type();
								if ( option[ backrub::detailed_balance ] ) {
									proposal_density_ratio = sidechainmover.last_proposal_density_ratio();
								}
								//TR << "proposal density: " << proposal_density_ratio << std::endl;
							}
						}
					}
				} else if ( move_prob > sc_prob ) {
					smallmover.apply(*pose);
					move_type = smallmover.type();
				} else {
					sidechainmover.apply(*pose);
					move_type = sidechainmover.type();
					if ( option[ backrub::detailed_balance ] ) {
						proposal_density_ratio = sidechainmover.last_proposal_density_ratio();
					}
				}

				//mdhist_proposed.record(pose->residue(1).chi());

				mc.boltzmann(*pose, move_type, proposal_density_ratio);

				//mdhist.record(pose->residue(1).chi());

				if ( option[ backrub::trajectory ] ) trajectory.update_after_boltzmann(mc);
				if ( option[ backrub::test_detailed_balance ] ) dof_recorder.update_after_boltzmann(*pose);
			}

			mc.show_counters();

			//utility::io::ozstream histstream(output_tag + "_hist.txt");
			//histstream << mdhist << mdhist_proposed;
			//histstream.close();

			if ( option[ backrub::test_detailed_balance ] ) {
				utility::io::ozstream histstream(output_tag + "_hist.txt");
				histstream << dof_recorder;
				histstream.close();

				utility::io::ozstream msestream(output_tag + "_mse.txt");
				dof_recorder.write_mse_summary(msestream);
				msestream.close();
			}

			// repack II: LAST
			// this could also be done IN the loop, so that we get a minimized structure right after the monte carlo run
			//packrotamersmover.apply(*pose);

			// dump out the low score and last accepted poses
			TR << "Last Score:" << std::endl;
			score_fxn->show(TR, *pose);
			TR.flush();

			*pose = mc.lowest_score_pose();

			// repack II: LOW
			//packrotamersmover.apply(*pose);

			TR << "Low Score:" << std::endl;
			score_fxn->show(TR, *pose);
			TR.flush();

			mc.last_accepted_pose().dump_pdb(output_tag + "_last.pdb");
			mc.lowest_score_pose().dump_pdb(output_tag + "_low.pdb");

			if ( custom_fold_tree ) {
				append_fold_tree_to_file(mc.last_accepted_pose().fold_tree(), output_tag + "_last.pdb");
				append_fold_tree_to_file(mc.lowest_score_pose().fold_tree(), output_tag + "_low.pdb");
			}
		}
	}

	// write parameters for any sets of branching atoms for which there were not optimization coefficients
	backrubmover.branchopt().write_database();

	return 0;
}
