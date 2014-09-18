// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

// Protocol Headers
#include <protocols/rigid/RigidBodyMover.hh>

// Core Headers
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



using namespace core;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, chain_for_destabilization )

static thread_local basic::Tracer TR( "apps.pilot.whitney_specificity_switch.main" );


/// General testing code
int
main( int argc, char * argv [] )
{

	try {


	NEW_OPT( chain_for_destabilization, "the chain ID which will harbor any mutations", "A" );

	devel::init(argc, argv);

	std::string const tmp_chain = option[ chain_for_destabilization ];
	if ( tmp_chain.length() != 1 ) {
		std::cerr << "ERROR!! Chain ID should be one character" << std::endl;
		exit(1);
	}
	char const chain_to_mutate = tmp_chain[0];

	TR << "Starting ddG calculations on chain " << chain_to_mutate << std::endl;

	// create pose for wild type native bound and unbound
	pose::Pose bound_pose_native, unbound_pose_native;

	//read in pdb file from command line
	std::string const input_pdb_name ( basic::options::start_file() );
	core::import_pose::pose_from_pdb( bound_pose_native, input_pdb_name );

	scoring::ScoreFunctionOP scorefxn( get_score_function() );
	scorefxn->set_weight( core::scoring::fa_dun, 0.1 );

	// Setup pose for unbound
	unbound_pose_native = bound_pose_native;

	core::Real const unbound_dist = 40.; // create unbound distance

	// rb_jump = rigid_backbone - changes the distance between the backbones of the molecules
  //   based upon the pdb file - so DON'T change numbers
	Size const rb_jump = 1;
	protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose_native, rb_jump );
	trans_mover.trans_axis( trans_mover.trans_axis() );
	trans_mover.step_size(unbound_dist);
	trans_mover.apply( unbound_pose_native );

	//write out pdb file for unbound native
	unbound_pose_native.dump_scored_pdb( "unbound_native.pdb", *scorefxn );

	// repacking unbound wt
	pack::task::PackerTaskOP base_packer_task_unbound( pack::task::TaskFactory::create_packer_task( unbound_pose_native) );
	base_packer_task_unbound->set_bump_check( false );
	base_packer_task_unbound->initialize_from_command_line();
	base_packer_task_unbound->or_include_current( true );


	// rpeack unbound native
	for (int ii = 1, nres = unbound_pose_native.total_residue(); ii < nres; ++ii ) {

		chemical::AA const native_aa ( unbound_pose_native.residue(ii).aa() );
		utility::vector1 <bool> repack_native( core::chemical::num_canonical_aas, false );
		repack_native[ native_aa ] = true;
		base_packer_task_unbound->nonconst_residue_task( ii ).restrict_to_repacking();

	}

	pack::pack_rotamers( unbound_pose_native, *scorefxn, base_packer_task_unbound );

	// calculate score for wild type native bound and unbound
	(*scorefxn)(bound_pose_native);
	(*scorefxn)(unbound_pose_native);

	// calculating wt bound and unbound for native
  core::Real const starting_bound_score = bound_pose_native.energies().total_energies()[ total_score ];
	core::Real const starting_unbound_score = unbound_pose_native.energies().total_energies()[ total_score ];

	TR << "Starting bound score for native is: " << starting_bound_score << std::endl;
	TR << "Starting unbound score for native is: " << starting_unbound_score << std::endl;

	// create delta G for native
	core::Real const dg_wt = starting_bound_score - starting_unbound_score;

	TR << "Delta G wild type score for native: " << dg_wt << std::endl;

	//set up copy of packer for bound native
  pack::task::PackerTaskOP base_packer_task_native( pack::task::TaskFactory::create_packer_task( bound_pose_native ) );
	base_packer_task_native->initialize_from_command_line();
	base_packer_task_native->set_bump_check( false );
	base_packer_task_native->or_include_current( true ); // allow crystallized position as possible rotamer

	// create and open file for printing scores
	utility::io::ozstream ddg_outstream;
	ddg_outstream.open( "ddg.out", std::ios::out );

	// create string for header for file
	std::ostringstream header_string_stream;
	header_string_stream << std::setw(9) << "Mutation";
	header_string_stream << std::setw(15) << "dg -wt";
	header_string_stream << std::setw(15) << "dg -mut";
	header_string_stream << std::setw(15) << "ddg";
	header_string_stream << std::setw(15) << "dg unbound";
	header_string_stream << std::setw(15) << "mut unbound";
	header_string_stream << std::setw(15) << "wt unbound";
	header_string_stream << std::setw(15) << "dg bound";
	header_string_stream << std::setw(15) << "mut bound";
	header_string_stream << std::setw(15) << "wt bound";
	ddg_outstream << header_string_stream.str() << std::endl;

	// Loop over all positions in native, sequentially mutate each residue to all possible
	// residues and calculate ddg in both the bound and the unbound
	for ( int j = 1, resnum = bound_pose_native.total_residue(); j <= resnum; ++j ) {

		// if the residue is not on desired chain, don't mutate
		if ( bound_pose_native.pdb_info()->chain(j) != chain_to_mutate ) continue;

		//set the packer to change aa at current position
		pack::task::PackerTaskOP position_packer_task_native( base_packer_task_native->clone() );
		utility::vector1 <bool>	allow_repacked_res( resnum, false );
		allow_repacked_res.at(j) = true;

		// find the neighbors and allow them to repack
		core::scoring::TenANeighborGraph const & graph = bound_pose_native.energies().tenA_neighbor_graph();
		for ( core::graph::Graph::EdgeListConstIter
						iter = graph.get_node( j )->const_edge_list_begin(),
						iter_end = graph.get_node( j )->const_edge_list_end();
					iter != iter_end; ++iter ) {
			Size const neighbor_id( (*iter)->get_other_ind( j ) );
			chemical::AA const neighbor_aa( bound_pose_native.residue(neighbor_id).aa() );
			if ( oneletter_code_from_aa(neighbor_aa) == 'C' ) continue;
			position_packer_task_native->nonconst_residue_task( neighbor_id ).restrict_to_repacking();
			allow_repacked_res.at(neighbor_id) = true;
			}

		position_packer_task_native->restrict_to_residues( allow_repacked_res );

	  // get the wild type amino acid at the current position
		chemical::AA const wt_aa( bound_pose_native.residue(j).aa() );

	  if ( oneletter_code_from_aa(wt_aa) == 'C' ) continue;

		for ( int res = 1; res <= core::chemical::num_canonical_aas; ++res ) {

			// if residue is a Cys, don't mutate
			if ( res == 2 ) continue;

			// create mutant poses, assign mutant poses to wild type
			pose::Pose mut_bound_native, mut_unbound_native;
			mut_bound_native = bound_pose_native;
			mut_unbound_native = unbound_pose_native;

			// Set up packer for bound, only allowing the current aa as mutant
			pack::task::PackerTaskOP repacked_packer_task_native( position_packer_task_native->clone() );
			utility::vector1 <bool> repack_aa( core::chemical::num_canonical_aas, false );
			repack_aa.at(res) = true;
			repacked_packer_task_native->nonconst_residue_task(j).restrict_absent_canonical_aas( repack_aa );

			// Mutate pose based upon selections
			pack::pack_rotamers( mut_bound_native, *scorefxn, repacked_packer_task_native );
			pack::pack_rotamers( mut_unbound_native, *scorefxn, repacked_packer_task_native );

			// calculate score after repacking
			(*scorefxn)(mut_bound_native);
			(*scorefxn)(mut_unbound_native);

			//determine current mutation - print to screen
			chemical::AA const mut_aa( mut_bound_native.residue(j).aa() );
			std::ostringstream mut_name;
		  mut_name << oneletter_code_from_aa(wt_aa) << " " << mut_bound_native.pdb_info()->number(j) << mut_bound_native.pdb_info()->chain(j) << " " << oneletter_code_from_aa( mut_aa );
			TR << "Processing Mutation: " << mut_name.str() << std::endl;

			// calculate new score after repacking and minimization
			(*scorefxn)(mut_bound_native);
			(*scorefxn)(mut_unbound_native);

			// calculate scores for mutant bound unbound
			core::Real const mut_bound_score = mut_bound_native.energies().total_energies()[ total_score ];
			core::Real const mut_unbound_score = mut_unbound_native.energies().total_energies()[ total_score ];

			// calculate ddg
			core::Real const dg_mut = mut_bound_score - mut_unbound_score;
			core::Real const ddg_native = dg_mut - dg_wt;

			// writes out new pdb file for mutant bound
			std::ostringstream fname;
			fname << oneletter_code_from_aa(wt_aa) << mut_bound_native.pdb_info()->number(j) << oneletter_code_from_aa( mut_aa );
			std::ostringstream outPDB_name;
			outPDB_name << fname.str() << ".pdb";
			mut_bound_native.dump_scored_pdb( outPDB_name.str(), *scorefxn );

			// Printing Results to screen
			TR << "Starting bound score for native: " << starting_bound_score << std::endl;
			TR << "Mutant bound score for native: " << mut_bound_score << std::endl;
			TR << "Delta G for mutant bound native: " <<  mut_bound_score - starting_bound_score << std::endl;

			TR << "Starting unbound score for native: " << starting_unbound_score << std::endl;
			TR << "Mutant unbound score for native: " << mut_unbound_score << std::endl;
			TR << "Delta G for mutant unbound: for native " << mut_unbound_score - starting_unbound_score << std::endl;


			// Reporting Data
			std::ostringstream data_string_stream;
			data_string_stream << std::setw(9) << mut_name.str();
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << dg_wt;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << dg_mut;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << ddg_native;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << mut_unbound_score - starting_unbound_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << mut_unbound_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << starting_unbound_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << mut_bound_score - starting_bound_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << mut_bound_score;
			data_string_stream << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(15) << starting_bound_score;

			// Print to file
			ddg_outstream << data_string_stream.str() << std::endl;

		}
	}

	// close file
	ddg_outstream.close();
	ddg_outstream.clear();

	TR << "Successfully finished computing ddGs" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



