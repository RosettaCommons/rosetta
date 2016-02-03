// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file cterdesign.cc
/// @brief backrub Monte Carlo/small phi/psi move/design protocol
/// @author Colin A. Smith (colin.smith@ucsf.edu)
/// @details

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/viewer/viewers.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/io/mpistream.hh>
#include <fstream>


static THREAD_LOCAL basic::Tracer TR( "apps.backrub" );

OPT_1GRP_KEY(Integer, backrub, ntrials)
OPT_1GRP_KEY(Real, backrub, sc_prob)
OPT_1GRP_KEY(Real, backrub, sm_prob)
OPT_1GRP_KEY(Real, backrub, sc_prob_uniform)
OPT_1GRP_KEY(Real, backrub, mc_kt)
OPT_1GRP_KEY(Real, backrub, mc_kt_initial)
OPT_1GRP_KEY(Real, backrub, mm_bend_weight)
OPT_1GRP_KEY(Integer, backrub, pack_frequency)
OPT_1GRP_KEY(StringVector, backrub, chains1)
OPT_1GRP_KEY(StringVector, backrub, chains2)

void *
my_main( void* );

int
main( int argc, char * argv [] )
{
    try {
	OPT(in::path::database);
	OPT(in::file::s);
	OPT(in::file::l);
	OPT(in::ignore_unrecognized_res);
	OPT(out::nstruct);
	OPT(packing::resfile);
	OPT(backrub::pivot_residues);
	OPT(backrub::pivot_atoms);
	OPT(backrub::min_atoms);
	OPT(backrub::max_atoms);
	NEW_OPT(backrub::ntrials, "number of Monte Carlo trials to run", 1000);
	NEW_OPT(backrub::sc_prob, "probability of making a side chain move", 0.25);
	NEW_OPT(backrub::sm_prob, "probability of making a small move", 0.25);
	NEW_OPT(backrub::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.1);
	NEW_OPT(backrub::mc_kt, "value of kT for Monte Carlo", 0.6);
	NEW_OPT(backrub::mc_kt_initial, "initial value of kT for Monte Carlo", 0.6);
	NEW_OPT(backrub::mm_bend_weight, "weight of mm_bend bond angle energy term", 1.0);
	NEW_OPT(backrub::pack_frequency, "number of moves between PackRotamers", 1000);
	NEW_OPT(backrub::chains1, "chains on one side of the interface", utility::vector1<std::string>(1, "A"));
	NEW_OPT(backrub::chains2, "chains on the other side of the interface", utility::vector1<std::string>(1, "B"));

	// initialize Rosetta
	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
	return 0;
}

core::scoring::TwoBodyEMapVector
interface_energymap(
	core::pose::Pose const & pose,
	std::set<core::Size> const & res_set1,
	std::set<core::Size> const & res_set2
)
{
	core::scoring::TwoBodyEMapVector energymap;
	core::scoring::EnergyGraph const & energygraph(pose.energies().energy_graph());

	// iterate over all residues in the first set
	for (std::set<core::Size>::const_iterator res_iter1 = res_set1.begin(); res_iter1 != res_set1.end(); ++res_iter1) {

		// get the node for the particular residue
		core::scoring::EnergyNode const * const node1(energygraph.get_energy_node(*res_iter1));

		// iterate over all edges connected to that node
		for (core::graph::EdgeListConstIterator edge_iter = node1->const_edge_list_begin();
		     edge_iter != node1->const_edge_list_end(); ++edge_iter) {

			// check to see if the other node connected to the edge is in the second set of residues
			if (res_set2.count((*edge_iter)->get_other_ind(*res_iter1))) {
				core::scoring::EnergyEdge const * const energyedge(dynamic_cast<core::scoring::EnergyEdge const *>(*edge_iter));
				assert(energyedge);
				// accumulate the energes in the overall TwoBodyEMapVector
				energymap += energyedge->fill_energy_map();
				//TR << "Adding scores for " << energyedge->get_first_node_ind() << "-" << energyedge->get_second_node_ind() << std::endl;
				//TR << energyedge->energy_map().show_nonzero() << std::endl;
			}
		}
	}

	return energymap;
}

core::scoring::TwoBodyEMapVector
interface_energymap(
	core::pose::Pose const & pose,
	std::set<char> const & chain_set1,
	std::set<char> const & chain_set2
)
{
	std::set<char> chains_intersect;
	std::set_intersection(chain_set1.begin(), chain_set1.end(), chain_set2.begin(), chain_set2.end(),
	                      std::insert_iterator<std::set<char> >(chains_intersect, chains_intersect.begin()));
	assert(chains_intersect.empty() );//size() == 0);

	std::set<core::Size> res_set1;
	std::set<core::Size> res_set2;

	for (core::Size i = 1; i <= pose.pdb_info()->nres(); ++i) {

		char chain = pose.pdb_info()->chain(i);
		if (chain_set1.count(chain)) {
			res_set1.insert(i);
		} else if (chain_set2.count(chain)) {
			res_set2.insert(i);
		}
	}
	/*
	TR << "Chain 1 Residues:";
	for (std::set<core::Size>::iterator iter = res_set1.begin(); iter != res_set1.end(); ++iter) {
		TR << " " << (*iter);
	}
	TR << std::endl;
	TR << "Chain 2 Residues:";
	for (std::set<core::Size>::iterator iter = res_set2.begin(); iter != res_set2.end(); ++iter) {
		TR << " " << (*iter);
	}
	TR << std::endl;
	*/

	return interface_energymap(pose, res_set1, res_set2);
}

void *
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// create a TaskFactory with the resfile
	using namespace core::pack::task;
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	}
	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( new operation::PreserveCBeta );

	// set up the score function and add the bond angle energy term
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
	score_fxn->set_weight(core::scoring::mm_bend, option[ backrub::mm_bend_weight ]);
	core::scoring::methods::EnergyMethodOptions energymethodoptions(score_fxn->energy_method_options());
	energymethodoptions.hbond_options()->decompose_bb_hb_into_pair_energies(true);
	energymethodoptions.bond_angle_central_atoms_to_score(option[ backrub::pivot_atoms ]);
	score_fxn->set_energy_method_options(energymethodoptions);

	// set up the BackrubMover
	protocols::backrub::BackrubMover backrubmover;
	// read known and unknown optimization parameters from the database
	backrubmover.branchopt().read_database();

	// set up the SmallMover
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->init_from_file(option[ in::file::movemap ]);
	protocols::simple_moves::SmallMover smallmover;
	smallmover.nmoves(1);
	smallmover.movemap(movemap);

	// set up the SidechainMover
	protocols::simple_moves::sidechain_moves::SidechainMover sidechainmover;
	sidechainmover.set_task_factory(main_task_factory);
	sidechainmover.set_prob_uniform(option[ backrub::sc_prob_uniform ]);

	// set up the PackRotamersMoveer
	protocols::simple_moves::PackRotamersMover packrotamersmover;
	packrotamersmover.task_factory(main_task_factory);
	packrotamersmover.score_function(score_fxn);

	utility::vector1<std::string> const & chain_vector1(option[ backrub::chains1 ]);
	utility::vector1<std::string> const & chain_vector2(option[ backrub::chains2 ]);
	std::set<char> chain_set1;
	std::set<char> chain_set2;
	for (utility::vector1<std::string>::const_iterator iter = chain_vector1.begin(); iter != chain_vector1.end(); ++iter) {
		if (iter->size() == 1) chain_set1.insert((*iter)[0]);
	}
	for (utility::vector1<std::string>::const_iterator iter = chain_vector2.begin(); iter != chain_vector2.end(); ++iter) {
		if (iter->size() == 1) chain_set2.insert((*iter)[0]);
	}

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	for (core::Size jobnum = 1; jobnum <= input_jobs.size(); ++jobnum) {

		TR << "Processing " << input_jobs[jobnum]->input_tag() << "..." << std::endl;

		// load the PDB file
		core::pose::PoseOP input_pose(new core::pose::Pose());
		if ( option[ in::file::centroid_input ].user() ) {
			TR.Warning << "*** This is untested with centroid mode! ***" << std::endl;
			core::import_pose::centroid_pose_from_pdb( *input_pose, input_jobs[jobnum]->input_tag() , core::import_pose::PDB_file);
		} else {
			core::import_pose::pose_from_file( *input_pose, input_jobs[jobnum]->input_tag() , core::import_pose::PDB_file);
		}

		backrubmover.clear_segments();
		backrubmover.set_input_pose(input_pose);

		TR << "Backrub segment lengths: " << option[ backrub::min_atoms ] << "-" << option[ backrub::max_atoms ] << " atoms"
		   << std::endl;

		TR << "Backrub main chain pivot atoms: " << option[ backrub::pivot_atoms ].value_string() << std::endl;

		backrubmover.add_mainchain_segments_from_options();

		TR << "Backrub Segments Added: " << backrubmover.num_segments() << std::endl;

		TR << "Score After PDB Load:" << std::endl;
		score_fxn->show(TR, *input_pose);
		TR.flush();

		backrubmover.optimize_branch_angles(*input_pose);
		sidechainmover.idealize_sidechains(*input_pose);

		TR << "Score After Branch Angle Optimization/Side Chain Idealization:" << std::endl;
		score_fxn->show(TR, *input_pose);
		TR.flush();

		protocols::moves::MonteCarlo mc(*input_pose, *score_fxn, option[ backrub::mc_kt ]);

		// create viewer windows if OpenGL is enabled
		protocols::viewer::add_monte_carlo_viewer(mc, "Backrub", 600, 600);

		// iterate to generate multiple structures
		for (int structnum = 1; structnum <= input_jobs[jobnum]->nstruct(); ++structnum)
		{
			// reset to the starting optimized pose
			core::pose::PoseOP pose(new core::pose::Pose(*input_pose));
			//mc.reset(*pose);

			// force a repack at the beginning
			packrotamersmover.apply(*pose);
			mc.reset(*pose);

			TR << "Score After Initial Pack:" << std::endl;
			score_fxn->show(TR, *pose);
			TR.flush();

			TR << "Running " << option[ backrub::ntrials ] << " trials..." << std::endl;

			for (int i = 1; i <= option[ backrub::ntrials ]; ++i) {

				std::string move_type;

				// could use random mover for this...
				core::Real move_prob = numeric::random::rg().uniform();
				if (i % option[ backrub::pack_frequency ] == 0) {
					packrotamersmover.apply(*pose);
					move_type = packrotamersmover.type();
				} else if (move_prob > option[ backrub::sm_prob ] + option[ backrub::sc_prob ]) {
					backrubmover.apply(*pose);
					move_type = backrubmover.type();
				} else if (move_prob > option[ backrub::sc_prob ]) {
					smallmover.apply(*pose);
					move_type = smallmover.type();
				} else {
					sidechainmover.apply(*pose);
					move_type = sidechainmover.type();
				}

				mc.boltzmann(*pose, move_type);
			}

			mc.show_counters();

			// dump out the low score and last accepted poses
			TR << "Last Score:" << std::endl;
			score_fxn->show(TR, *pose);
			TR.flush();
			core::scoring::TwoBodyEMapVector last_interface_energies(interface_energymap(*pose, chain_set1, chain_set2));
			TR << "Last Interface Components:";
			last_interface_energies.show_weighted(TR, score_fxn->weights());
			TR << std::endl;
			TR << "Last Interface Total: " << last_interface_energies.dot(score_fxn->weights()) << std::endl;

			std::string output_tag(input_jobs[jobnum]->output_tag(structnum));

			pose->dump_pdb(output_tag + "_last.pdb");
			std::ofstream last_outfile((output_tag + "_last.pdb").c_str(), std::ios::out|std::ios::app);
			score_fxn->show(last_outfile, *pose);
			last_outfile << "Interface Components:";
			last_interface_energies.show_weighted(last_outfile, score_fxn->weights());
			last_outfile << std::endl;
			last_outfile << "Interface Total: " << last_interface_energies.dot(score_fxn->weights()) << std::endl;
			last_outfile.close();

			*pose = mc.lowest_score_pose();

			TR << "Low Score:" << std::endl;
			score_fxn->show(TR, *pose);
			TR.flush();
			core::scoring::TwoBodyEMapVector low_interface_energies(interface_energymap(*pose, chain_set1, chain_set2));
			TR << "Low Interface Components:";
			low_interface_energies.show_weighted(TR, score_fxn->weights());
			TR << std::endl;
			TR << "Low Interface Total: " << low_interface_energies.dot(score_fxn->weights()) << std::endl;

			pose->dump_pdb(output_tag + "_low.pdb");
			std::ofstream low_outfile((output_tag + "_low.pdb").c_str(), std::ios::out|std::ios::app);
			score_fxn->show(low_outfile, *pose);
			low_outfile << "Interface Components:";
			low_interface_energies.show_weighted(low_outfile, score_fxn->weights());
			low_outfile << std::endl;
			low_outfile << "Interface Total: " << low_interface_energies.dot(score_fxn->weights()) << std::endl;
			low_outfile.close();
		}
	}

	// write parameters for any sets of branching atoms for which there were not optimization coefficients
	backrubmover.branchopt().write_database();

	return 0;
}
