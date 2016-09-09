// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <protocols/rigid/RigidBodyMover.hh>

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//OPT_KEY( String, chain_for_destabilization )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.nadeem_find_destabilizing_mutation.main" );

int
main( int argc, char * argv [] )
{
	try {

		//NEW_OPT( chain_for_destabilization, "the chain ID which will harbor the destabilizing mutation", "A" );

		devel::init(argc, argv);
		/*
		std::string const tmp_chain = option[ chain_for_destabilization ];
		if ( tmp_chain.length() != 1 ) {
		std::cerr << "ERROR!! Chain ID should be one character" << std::endl;
		exit(1);
		}
		char const chain_to_analyze = tmp_chain[0];
		*/

		// Number of times to repack the starting structure (to look for global minimum)
		//mjo commenting out 'num_starting_repacks' because it is unused and causes a warning
		//Size const num_starting_repacks(500);
		// Number of times to repack each mutant
		Size const num_repacks(20);

		pose::Pose wt_pose;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( wt_pose, input_pdb_name , core::import_pose::PDB_file);

		// Setup for scoring/repacking
		scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(SOFT_REP_DESIGN_WTS) );
		// scorefxn->set_weight( core::scoring::fa_dun, 0.1 );
		(*scorefxn)(wt_pose);

		// Repack the WT
		pack::task::PackerTaskOP repack_task( pack::task::TaskFactory::create_packer_task( wt_pose ));
		repack_task->set_bump_check( false );
		repack_task->initialize_from_command_line();
		repack_task->or_include_current( true );
		for ( int resnum=1, totres = wt_pose.size(); resnum <= totres; ++resnum ) {
			repack_task->nonconst_residue_task( resnum ).restrict_to_repacking();
		}

		pack::pack_rotamers_loop( wt_pose, *scorefxn, repack_task, 10 * num_repacks );
		(*scorefxn)(wt_pose);
		core::Real const wt_score = wt_pose.energies().total_energies()[ total_score ];

		// For each jump, setup the unbound pose and compute dG
		/* utility::vector1< pose::Pose > unbound_poses( wt_pose.num_jump() );
		utility::vector1< core::Real > wt_dG_binding( wt_pose.num_jump() );
		for ( core::Size j=1; j <= wt_pose.num_jump(); ++j) {
		unbound_poses[j] = wt_pose;
		core::Real const unbound_dist = 40.;
		protocols::rigid::RigidBodyTransMover wt_trans_mover( wt_pose, j );
		wt_trans_mover.trans_axis( wt_trans_mover.trans_axis() );
		wt_trans_mover.step_size(unbound_dist);
		wt_trans_mover.apply( unbound_poses[j] );
		(*scorefxn)(unbound_poses[j]);
		pack::pack_rotamers_loop( unbound_poses[j], *scorefxn, repack_task, 10 * num_repacks );
		(*scorefxn)(unbound_poses[j]);
		wt_dG_binding[j] = wt_bound_score - (unbound_poses[j]).energies().total_energies()[ total_score ];
		(unbound_poses[j]).dump_scored_pdb( "unbound"+utility::to_string(j)+".pdb", *scorefxn );
		}
		*/
		// Open output file, generate the header line (save it for printing in the log later), print to file
		std::stringstream filename;
		if ( !option[ OptionKeys::out::output_tag ]().empty() ) {
			filename<<option[ OptionKeys::out::output_tag ]()<<".mutations.out";
		} else {
			filename<<basic::options::start_file()<<".mutations.out";
		}

		utility::io::ozstream outstream;
		outstream.open(filename.str(), std::ios::out);
		//outstream << "Mutation poly_ala_rep ddG_jump1 ddG_jump2 ...." << std::endl;

		//mjo commenting out 'aa_cys' because it is unused and causes a warning
		//int const aa_cys = 2;

		//mjo commenting out 'aa_gly' because it is unused and causes a warning
		//int const aa_gly = 6;

		//mjo commenting out 'aa_pro' because it is unused and causes a warning
		//int const aa_pro = 13;

		//mjo commenting out 'ala_polyA_farep' because it is unused and causes a warning
		//core::Real ala_polyA_farep = -999.;

		// loop over all residues on the desired chain, all amino acids except C, G, P
		for ( int resnum = 1,totres =  wt_pose.size(); resnum <= totres; ++resnum ) {
			//if ( wt_pose.pdb_info()->chain(resnum) != chain_to_analyze ) continue;
			for ( int aa = 1; aa <= core::chemical::num_canonical_aas; ++aa ) {
				/*
				if ( aa == aa_cys ) continue;
				if ( aa == aa_gly ) continue;
				if ( aa == aa_pro ) continue;
				*/
				chemical::AA const wt_aa( wt_pose.residue(resnum).aa());
				/*   if ( oneletter_code_from_aa(wt_aa) == 'C' ) continue;
				if ( oneletter_code_from_aa(wt_aa) == 'G' ) continue;
				if ( oneletter_code_from_aa(wt_aa) == 'P' ) continue;
				*/
				// JK DEBUG
				if ( aa == wt_aa ) continue;

				// create poses for mutants
				pose::Pose mut_bound;
				mut_bound = wt_pose;
				//   pose::Pose mut_polyA;
				//   mut_polyA = wt_pose;
				(*scorefxn)(mut_bound);
				//   (*scorefxn)(mut_polyA);

				pack::task::PackerTaskOP mut_task( pack::task::TaskFactory::create_packer_task( mut_bound ));
				mut_task->set_bump_check( false );
				mut_task->initialize_from_command_line();
				mut_task->or_include_current( true );
				/*
				pack::task::PackerTaskOP polyA_packer_task( pack::task::TaskFactory::create_packer_task( mut_polyA ));
				polyA_packer_task->set_bump_check( false );
				polyA_packer_task->initialize_from_command_line();
				polyA_packer_task->or_include_current( true );
				*/
				// restrict packer task to single sequence position of interest
				utility::vector1<bool> allow_repacked( totres, false );
				allow_repacked.at(resnum) = true;
				utility::vector1< bool > ala_only( core::chemical::num_canonical_aas, false );
				ala_only.at(1) = true;

				// code for repacking neighbors...
				core::scoring::TenANeighborGraph const & graph = mut_bound.energies().tenA_neighbor_graph();
				for ( utility::graph::Graph::EdgeListConstIter
						iter = graph.get_node( resnum )->const_edge_list_begin(),
						iter_end = graph.get_node( resnum )->const_edge_list_end();
						iter != iter_end; ++iter ) {
					Size const neighbor_res( (*iter)->get_other_ind( resnum ) );
					mut_task->nonconst_residue_task( neighbor_res ).restrict_to_repacking();
					allow_repacked.at(neighbor_res) = true;
					//    polyA_packer_task->nonconst_residue_task( neighbor_res ).restrict_absent_canonical_aas( ala_only );
				}
				mut_task->restrict_to_residues( allow_repacked );
				//   polyA_packer_task->restrict_to_residues( allow_repacked );

				// set residue to allow in packer task
				utility::vector1< bool > mut_resnumlist( core::chemical::num_canonical_aas, false );
				mut_resnumlist.at(aa) = true;
				mut_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas( mut_resnumlist );
				//    polyA_packer_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas( mut_resnumlist );

				// make the mutation to mutant bound
				pack::pack_rotamers_loop( mut_bound, *scorefxn, mut_task, num_repacks );
				(*scorefxn)(mut_bound);
				core::Real const mutant_score = mut_bound.energies().total_energies()[ total_score ];

				std::ostringstream mutation_name;
				chemical::AA const mut_aa( mut_bound.residue(resnum).aa());
				mutation_name << oneletter_code_from_aa(wt_aa) << mut_bound.pdb_info()->number(resnum) << oneletter_code_from_aa(mut_aa);
				//   mut_bound.dump_scored_pdb("mutant_"+mutation_name.str()+".pdb", *scorefxn);

				// incorporate this mutation into a poly-Ala scaffold, find how much MORE it clashes than Ala does
				//   pack::pack_rotamers_loop( mut_polyA, *scorefxn, polyA_packer_task, num_repacks );
				//   (*scorefxn)(mut_polyA);
				//   if ( aa == 1 ) {
				//    ala_polyA_farep = mut_polyA.energies().total_energies()[ fa_rep ];
				//   }
				//   core::Real mut_polyA_farep = mut_polyA.energies().total_energies()[ fa_rep ] - ala_polyA_farep;

				// Print to file
				outstream << mutation_name.str() << " " << mutant_score - wt_score<<std::endl;

				/*   for ( core::Size j=1; j <= wt_pose.num_jump(); ++j) {
				pose::Pose mut_unbound = unbound_poses[j];
				(*scorefxn)(mut_unbound);
				pack::pack_rotamers_loop( mut_unbound, *scorefxn, mut_task, num_repacks );
				(*scorefxn)(mut_unbound);
				core::Real const mutant_unbound_score = mut_unbound.energies().total_energies()[ total_score ];
				core::Real const dGbinding_mut = mutant_score - mutant_unbound_score;
				core::Real const ddG = dGbinding_mut - wt_dG_binding[j];
				outstream <<" "<< ddG;
				*/
				// JK DEBUG
				//    outstream <<" "<< mutant_score - wt_score;
				//    outstream <<" "<< mutant_unbound_score - (unbound_poses[j]).energies().total_energies()[ total_score ];

				//   }

				//   outstream << std::endl;
			}
		}

		outstream.close();
		outstream.clear();

		TR << "Successfully finished computing ddGs" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


