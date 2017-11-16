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
#include <basic/options/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "apps.pilot.nadeem_compensatory_mutation.main" );


int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		TR << "Starting ddG calculations" << std::endl;

		//Take wild type NGF from pdb
		TR << "About to read bound pose" << std::endl;
		pose::Pose wt_NGFbound;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( wt_NGFbound, input_pdb_name , core::import_pose::PDB_file);

		// Score for the wildtype
		// Setup for scoring/repacking
		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		scorefxn->set_weight( core::scoring::fa_dun, 0.1 );

		(*scorefxn)(wt_NGFbound);

		// Open output file, generate the header line (save it for printing in the log later), print to file
		std::string ddg_outfname = "compensatory_mutations.out";
		utility::io::ozstream ddg_outstream;
		ddg_outstream.open(ddg_outfname, std::ios::out);
		ddg_outstream << "Base_mutation total_score delta_score" << std::endl;

		int const aa_cys = 2;
		int const aa_gly = 6;
		int const aa_pro = 13;

		// loop over all residues, all amino acids
		for ( int resnum = 1,totres =  wt_NGFbound.size(); resnum <= totres; ++resnum ) {

			// NA DEBUG
			//  if ( ( wt_NGFbound.pdb_info()->number(resnum) != 54 ) ) continue;
			//  if ( ( wt_NGFbound.pdb_info()->chain(resnum) != 'V' ) ) continue;

			for ( int aa = 1; aa <= core::chemical::num_canonical_aas; ++aa ) {

				if ( aa == aa_cys ) continue;
				if ( aa == aa_gly ) continue;
				if ( aa == aa_pro ) continue;

				chemical::AA const wt_aa( wt_NGFbound.residue(resnum).aa());
				if ( oneletter_code_from_aa(wt_aa) == 'C' ) continue;
				if ( oneletter_code_from_aa(wt_aa) == 'G' ) continue;
				if ( oneletter_code_from_aa(wt_aa) == 'P' ) continue;

				// create poses for mutant
				pose::Pose mut_NGFbound;
				mut_NGFbound = wt_NGFbound;
				(*scorefxn)(mut_NGFbound);

				pack::task::PackerTaskOP mutation_packer_task( pack::task::TaskFactory::create_packer_task( mut_NGFbound ));
				mutation_packer_task->set_bump_check( false );
				mutation_packer_task->initialize_from_command_line();
				mutation_packer_task->or_include_current( true );

				// restrict packer task to single sequence position of interest
				utility::vector1<bool> allow_repacked( totres, false );
				allow_repacked.at(resnum) = true;

				// residues to allow in redesign - all but CGP
				utility::vector1< bool > no_CGP_redesign( core::chemical::num_canonical_aas, true );
				no_CGP_redesign.at(aa_cys) = false;
				no_CGP_redesign.at(aa_gly) = false;
				no_CGP_redesign.at(aa_pro) = false;

				// code for repacking neighbors...
				core::scoring::TenANeighborGraph const & graph = mut_NGFbound.energies().tenA_neighbor_graph();
				for ( utility::graph::Graph::EdgeListConstIter
						iter = graph.get_node( resnum )->const_edge_list_begin(),
						iter_end = graph.get_node( resnum )->const_edge_list_end();
						iter != iter_end; ++iter ) {
					Size const neighbor_res( (*iter)->get_other_ind( resnum ) );
					mutation_packer_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas( no_CGP_redesign );
					allow_repacked.at(neighbor_res) = true;
				}
				mutation_packer_task->restrict_to_residues( allow_repacked );

				// set residue to allow in packer task
				utility::vector1< bool > repack_resnumlist( core::chemical::num_canonical_aas, false );
				repack_resnumlist.at(aa) = true;
				mutation_packer_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas( repack_resnumlist );

				// make the mutation to the bound
				pack::pack_rotamers( mut_NGFbound, *scorefxn, mutation_packer_task );

				// find the score for the bound
				(*scorefxn)(mut_NGFbound);
				core::Real const mutation_bound_score = mut_NGFbound.energies().total_energies()[ total_score ];
				TR << "MUT bound score is: " << mutation_bound_score << std::endl;

				std::ostringstream base_mutation_name;
				chemical::AA const mut_aa( mut_NGFbound.residue(resnum).aa());
				base_mutation_name << oneletter_code_from_aa(wt_aa) << mut_NGFbound.pdb_info()->number(resnum) << oneletter_code_from_aa(mut_aa) << "_chain" << mut_NGFbound.pdb_info()->chain(resnum);
				TR << "Base mutation: " << base_mutation_name.str() << std::endl;
				mut_NGFbound.dump_scored_pdb("OUTPUT_PDBs/redesign_"+base_mutation_name.str()+".pdb", *scorefxn);

				// loop over all residues
				for ( int mutres = 1; mutres <= totres; ++mutres ) {
					if ( mutres == resnum ) continue;
					chemical::AA const orig_aa( wt_NGFbound.residue(mutres).aa());
					chemical::AA const mut_aa( mut_NGFbound.residue(mutres).aa());
					if ( mut_aa != orig_aa ) {
						std::ostringstream compensatory_mutation_name;
						compensatory_mutation_name << oneletter_code_from_aa(orig_aa) << mut_NGFbound.pdb_info()->number(mutres) << oneletter_code_from_aa(mut_aa) << "_chain" << mut_NGFbound.pdb_info()->chain(mutres);
						TR << "Compensatory mutation: " << compensatory_mutation_name.str() << std::endl;
					}
				}

				// Build the unbound, by pulling apart the redesign
				TR << "About to build the unbound" << std::endl;
				pose::Pose mut_NGFunbound;
				mut_NGFunbound = mut_NGFbound;
				core::Real const unbound_dist = 40.;
				Size const rb_jump = 1; // use the first jump as the one between partners
				protocols::rigid::RigidBodyTransMover wt_trans_mover( mut_NGFunbound, rb_jump );
				wt_trans_mover.trans_axis( wt_trans_mover.trans_axis() );
				wt_trans_mover.step_size(unbound_dist);
				wt_trans_mover.apply( mut_NGFunbound );
				(*scorefxn)(mut_NGFunbound);

				// Repack the unbound
				TR << "About to repack the unbound" << std::endl;
				pack::task::PackerTaskOP unbound_packer_task( pack::task::TaskFactory::create_packer_task( mut_NGFunbound ));
				unbound_packer_task->set_bump_check( false );
				unbound_packer_task->initialize_from_command_line();
				unbound_packer_task->or_include_current( true );
				for ( int resnum =1, totres = mut_NGFunbound.size(); resnum <= totres; ++resnum ) {
					unbound_packer_task->nonconst_residue_task( resnum ).restrict_to_repacking();
				}
				unbound_packer_task->restrict_to_residues( allow_repacked );
				pack::pack_rotamers( mut_NGFunbound, *scorefxn, unbound_packer_task );
				(*scorefxn)(mut_NGFunbound);

				// find the score differences for the mutant unbound
				core::Real const mutation_unbound_score = mut_NGFunbound.energies().total_energies()[ total_score ];
				TR << "MUT unbound score is: " << mutation_unbound_score << std::endl;
				// find the mutant delta_score
				core::Real const mut_delta_score = mutation_bound_score - mutation_unbound_score;
				TR << "MUT delta_score is: " << mut_delta_score << std::endl;

				ddg_outstream << base_mutation_name.str() << " " << std::setiosflags(std::ios::fixed) << mutation_bound_score<< " " << std::setiosflags(std::ios::fixed) << mut_delta_score << std::endl;
			}
		}

		ddg_outstream.close();
		ddg_outstream.clear();

		TR << "Successfully finished redesigning compensatory mutations" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

