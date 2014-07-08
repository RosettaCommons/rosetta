// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ig_dump.cc
/// @brief tool to calculate and dump out an interaction graph
/// @author Colin A. Smith (colin.smith@ucsf.edu)
/// @detailed
/// For every PDB file given as input, a directory is created containing the interaction graph data.
/// It contains two types of tab-delimited, gzipped files:
///
/// <RESNUM>.gz: These files contain the one body energies for each rotamer at the given residue number.
/// The colums are as follows:
///
/// 1: one body energy
/// 2: amino acid type
/// 3-7: values of the sidechain dihedral angles 1-4
///
/// <RESNUM1>-<RESNUM2>.gz: These files contain the two body energy matricies for rotamers at the given
/// residue numbers. There is a row for every rotamer at RESNUM1 and a column for every rotamer at
/// RESNUM2.

// Protocols Headers
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// Platform Headers
#include <platform/types.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>




basic::Tracer TR("apps.ig_dump");

int
main( int argc, char * argv [] )
{
    try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// initialize Rosetta
	devel::init(argc, argv);

	// create a TaskFactory with the resfile
	using namespace core::pack::task;
	TaskFactoryOP main_task_factory = new TaskFactory;
	main_task_factory->push_back( new operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new operation::ReadResfile );
	}

	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
	core::pose::PoseOP input_pose;

	for (core::Size i = 1; i <= input_jobs.size(); ++i) {

		TR << "Processing " << input_jobs[i]->input_tag() << "..." << std::endl;

		// load the PDB file
		input_pose = new core::pose::Pose();
		if ( option[ in::file::centroid_input ].user() ) {
			core::import_pose::centroid_pose_from_pdb( *input_pose, input_jobs[i]->input_tag() );
		} else {
			core::import_pose::pose_from_pdb( *input_pose, input_jobs[i]->input_tag() );
		}

		std::string output_tag(input_jobs[i]->output_tag(1));

		// make a directory for the results if necessary
		if ( utility::file::file_exists(output_tag) ) {
			TR << "Using existing " << output_tag << " directory for interaction graph data" << std::endl;
		} else {
			TR << "Creating " << output_tag << " directory for interaction graph data" << std::endl;
			if (  !utility::file::create_directory(output_tag) ) {
				utility_exit_with_message("Could not create results directory");
			}
		}

		// allocate variables necessary for creating an interaction graph
		PackerTaskCOP task = main_task_factory->create_task_and_apply_taskoperations( *input_pose );
		core::pack::rotamer_set::RotamerSetsOP rotsets = new core::pack::rotamer_set::RotamerSets;
		core::pack::interaction_graph::InteractionGraphBaseOP ig;

		// create rotamers and calculate interaction graph energies
		core::pack::pack_rotamers_setup(*input_pose, *score_fxn, task, rotsets, ig);

		int const num_nodes(ig->get_num_nodes());

		TR << "One Body Energies:" << std::endl << std::endl;

		for (int i = 1; i <= num_nodes; ++i) {

			TR << "Writing Node " << rotsets->moltenres_2_resid(i) << std::endl;

			std::ostringstream fname;
			fname << output_tag << platform::file::PATH_SEPARATOR << rotsets->moltenres_2_resid(i) << ".gz";

			utility::io::ozstream outfile(fname.str());

			int const num_states(ig->get_num_states_for_node(i));
			for (int j = 1; j <= num_states; ++j) {
				core::conformation::Residue const & rotres(*rotsets->rotamer_set_for_moltenresidue(i)->rotamer(j));
				outfile << ig->get_one_body_energy_for_node_state(i, j) << "\t" << rotres.name3() ;
				for (core::Size k = 1; k <= rotres.nchi(); ++k) {
					outfile << "\t" << rotres.chi(k);
				}
				outfile << "\n"; // don't flush with std::endl!
			}

			outfile.close();
		}
		TR << std::endl;
		TR << "Two Body Energies:" << std::endl << std::endl;

		core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig =
		dynamic_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph * > ( ig.get() );
		if ( !pig ) {
			utility_exit_with_message("Interaction graph is not pre-computed");
		}

		for (int node1 = 1; node1 <= num_nodes; ++node1) {
			for (int node2 = node1+1; node2 <= num_nodes; ++node2) {

				int const num_states1(ig->get_num_states_for_node(node1));
				int const num_states2(ig->get_num_states_for_node(node2));

				if (!ig->get_edge_exists(node1, node2)) continue;

				TR << "Writing Node " << rotsets->moltenres_2_resid(node1) << " - Node " << rotsets->moltenres_2_resid(node2)
				   << std::endl;

				std::ostringstream fname;
				fname << output_tag << platform::file::PATH_SEPARATOR << rotsets->moltenres_2_resid(node1) << "-"
				      << rotsets->moltenres_2_resid(node2) << ".gz";

				utility::io::ozstream outfile(fname.str());

				for (int state1 = 1; state1 <= num_states1; ++state1) {
					for (int state2 = 1; state2 <= num_states2; ++state2) {

						if (state2 > 1) outfile << "\t";
						outfile << pig->get_two_body_energy_for_edge(node1, node2, state1, state2);
					}
					outfile << "\n"; // don't flush with std::endl!
				}

				outfile.close();
			}
		}
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
       return 0;
}
