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

#include <devel/init.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, phosphotyr_num )
OPT_KEY( String, phosphotyr_chain )
OPT_KEY( Real, match_distance_cutoff )
OPT_KEY( Real, phosphate_force_constant )
OPT_KEY( Boolean, do_minimization )

static basic::Tracer TR( "apps.pilot.phosphoArg_build_rotamers.main" );


int
main( int argc, char * argv [] )
{

	NEW_OPT( phosphotyr_num, "which residue the starting pY is", 0 );
	NEW_OPT( phosphotyr_chain, "which chain is the starting pY is on", "P" );
	NEW_OPT( match_distance_cutoff, "required distance from Arg N to pTyr O", 1.0 );
	NEW_OPT( phosphate_force_constant, "force constant for the coordinate constraint used in minimization", 25.0 );
	NEW_OPT( do_minimization, "whether or not to do minimzation for rotamers found", true );

	devel::init(argc, argv);

	TR << "Starting phospho-Arg calculations" << std::endl;

	std::string const tmp_chain = option[ phosphotyr_chain ];
	if ( tmp_chain.length() != 1 ) {
		std::cerr << "ERROR!! Chain ID should be one character" << std::endl;
		exit(1);
	}
	char const phosphotyr_chain_char = tmp_chain[0];

	core::Real const force_constant = option[ phosphate_force_constant ];
	core::Real const match_dist_cut = option[ match_distance_cutoff ];

	TR << "Reading pose" << std::endl;

	pose::Pose input_pose;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( input_pose, input_pdb_name );

	//	if ( input_pose.residue( input_pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
	//		input_pose.append_residue_by_jump
	//			( *conformation::ResidueFactory::create_residue( input_pose.residue(1).residue_type_set().name_map( "VRT" ) ),
	//				input_pose.total_residue()/2 );
	//	}

	Size const totres = input_pose.total_residue();
	pose::Pose pose = input_pose;

	// Setup for scoring/repacking
	scoring::ScoreFunctionOP scorefxn( scoring::getScoreFunction() );
	scorefxn->set_weight( core::scoring::fa_dun, 0. );
	scorefxn->set_weight( core::scoring::fa_intra_rep, 0. );
	scorefxn->set_weight( core::scoring::coordinate_constraint, 0.5 );
	(*scorefxn)(pose);

	// The resnum for the phosphoTyr which will be replaced by phosphoArg
	int const pY_pdb_number = option[ phosphotyr_num ];

	Size pR_resnum = 0;
 	for ( Size resnum = 1; resnum <= totres; ++resnum ) {
		if ( ( pose.pdb_info()->number(resnum) == pY_pdb_number ) &&
			( pose.pdb_info()->chain(resnum) == phosphotyr_chain_char ) ) {
			if ( pR_resnum != 0 ) {
				std::cerr << "ERROR!! Found more than one possible phosphoTyr" << std::endl;
				exit(1);
			}
			pR_resnum = resnum;
		}
	}
	if ( pR_resnum == 0 ) {
		std::cerr << "ERROR!! Could not find phosphoTyr" << std::endl;
		exit(2);
	}
	// store the location of critical atoms in the starting Tyr
	conformation::Atom const pTyr_P = pose.residue(pR_resnum).atom( "P" );
	Vector const pTyr_P_xyz = pTyr_P.xyz();
	conformation::Atom const pTyr_OH = pose.residue(pR_resnum).atom( "OH" );
	conformation::Atom const pTyr_O1P = pose.residue(pR_resnum).atom( "O1P" );
	conformation::Atom const pTyr_O2P = pose.residue(pR_resnum).atom( "O2P" );
	conformation::Atom const pTyr_O3P = pose.residue(pR_resnum).atom( "O3P" );

	// build a poly-Ala pose
	pack::task::PackerTaskOP polyA_task( pack::task::TaskFactory::create_packer_task( pose ));
	polyA_task->set_bump_check( false );
	polyA_task->initialize_from_command_line();
	polyA_task->or_include_current( false );

	int const aa_ala = 1;
	utility::vector1< bool > redesign_ala( core::chemical::num_canonical_aas, false );
	redesign_ala.at(aa_ala) = true;

	utility::vector1<bool> polyA_allow_repacked( totres, true );
	polyA_allow_repacked.at(pR_resnum) = false;
	polyA_task->restrict_to_residues( polyA_allow_repacked );

 	for ( Size resnum = 1; resnum <= totres; ++resnum ) {
		if ( resnum == pR_resnum ) continue;
		chemical::AA const aa( pose.residue(resnum).aa());
		if ( ( oneletter_code_from_aa(aa) == 'G' ) || ( oneletter_code_from_aa(aa) == 'P' ) ) {
			// preserve Gly and Pro
			polyA_task->nonconst_residue_task(resnum).restrict_to_repacking();
		} else {
			// change rest to Ala
			polyA_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas( redesign_ala );
		}
	}
	pack::pack_rotamers( pose, *scorefxn, polyA_task );
	(*scorefxn)(pose);

	// build rotamers for Arg instead of the phosphoTyr
	pack::task::PackerTaskOP arg_task( pack::task::TaskFactory::create_packer_task( pose ));
	arg_task->set_bump_check( true );
	arg_task->initialize_from_command_line();
	arg_task->or_include_current( false );

	// restrict packer task to single sequence position of interest
	utility::vector1<bool> allow_repacked( totres, false );
	allow_repacked.at(pR_resnum) = true;
	arg_task->restrict_to_residues( allow_repacked );

	int const aa_arg = 15;
	utility::vector1< bool > redesign_arg( core::chemical::num_canonical_aas, false );
	redesign_arg.at(aa_arg) = true;
	arg_task->nonconst_residue_task(pR_resnum).restrict_absent_canonical_aas( redesign_arg );
	arg_task->nonconst_residue_task(pR_resnum).or_ex1( true );
	arg_task->nonconst_residue_task(pR_resnum).or_ex2( true );
	arg_task->nonconst_residue_task(pR_resnum).or_ex3( true );
	arg_task->nonconst_residue_task(pR_resnum).or_ex4( true );
	arg_task->nonconst_residue_task(pR_resnum).and_extrachi_cutoff( 0 );

	pack::pack_scorefxn_pose_handshake( pose, *scorefxn);
	pose.update_residue_neighbors();
	scorefxn->setup_for_packing( pose, arg_task->repacking_residues(), arg_task->designing_residues()  );

	//replace this with RotSetsFactory
	pack::rotamer_set::RotamerSetsOP rotsets( new pack::rotamer_set::RotamerSets() );
	rotsets->set_task( arg_task );
	graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, arg_task );
	rotsets->build_rotamers( pose, *scorefxn, packer_neighbor_graph );
	//	rotsets->dump_pdb( pose, basic::options::start_file()+".all_arg_rotamers.pdb" );

	// pull out the RotamerSet for the "interesting" position
	pack::rotamer_set::RotamerSetOP pArg_rotset = rotsets->rotamer_set_for_residue( pR_resnum );
	TR << "Built " << pArg_rotset->num_rotamers() << " Arg rotamers" << std::endl;

	// find those where the Arg has an outer N close to the Tyr O
	std::string const out_rotamer_fname( basic::options::start_file()+".good_arg_rotamers.pdb" );
	std::ofstream out_stream( out_rotamer_fname.c_str() );

	Size model(0), atom_counter(0);
 	for ( Size r = 1; r <= pArg_rotset->num_rotamers(); ++r ) {
		//		conformation::ResidueOP rotamer = pArg_rotset->nonconst_rotamer(r); // jk if needed....
		conformation::ResidueCOP rotamer = pArg_rotset->rotamer(r);

		if ( ! rotamer->has( "P" ) ) continue;
		conformation::Atom const pArg_P = rotamer->atom( "P" );
		if ( distance(pArg_P.xyz(),pTyr_P_xyz) > match_dist_cut ) continue;

		++model;
		out_stream << "MODEL     " << model << "\n";
		// JK DEBUG
		//		io::pdb::dump_pdb_residue( *rotamer, atom_counter, out_stream );
		out_stream << "ENDMDL\n";

		pose::Pose min_pose = pose;

		// put this rotamer onto input_pose (replacing Tyr)
		min_pose.replace_residue( pR_resnum, *rotamer, false );
		// JK DEBUG
		min_pose.dump_scored_pdb( basic::options::start_file()+".premin"+utility::to_string(model)+".pdb", *scorefxn );

		if ( option[ do_minimization ] ) {

			// restore the sidechains from input pose, except the one we minimized
			for ( Size resnum = 1; resnum <= totres; ++resnum ) {
				if ( resnum == pR_resnum ) continue;
				min_pose.replace_residue( resnum, input_pose.residue(resnum), false );
			}

			// carry out minimzation of (only) Arg chi angles
			kinematics::MoveMap mm_pArg;
			mm_pArg.set_chi( false );
			mm_pArg.set_bb( false );
			mm_pArg.set_jump( false );
			mm_pArg.set_chi( pR_resnum, true );

			// apply constraint to pose - P and N atoms only (for now)
			id::AtomID P_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("P") , pR_resnum );
			id::AtomID fixed_atom_id = id::AtomID( min_pose.residue(1).atom_index("CA") , 1 );
			//			id::AtomID fixed_atom_id = id::AtomID( 1 , min_pose.total_residue() );

			core::scoring::constraints::FuncOP func( new scoring::constraints::HarmonicFunc( 0., 1./force_constant ) );
			min_pose.add_constraint( new core::scoring::constraints::CoordinateConstraint( P_atom_id, fixed_atom_id, pTyr_P_xyz, func ) );

			// apply constraint to bring N of pArg onto OH of pTyr
			// figure out which pArg N is closest to the phosphate, rather than figure out connectivity
			conformation::Atom const pArg_P = min_pose.residue(pR_resnum).atom( "P" );
			core::Real const distN1 = distance(min_pose.residue(pR_resnum).atom( "NH1" ).xyz(),pArg_P.xyz() );
			core::Real const distN2 = distance(min_pose.residue(pR_resnum).atom( "NH2" ).xyz(),pArg_P.xyz() );
			bool N1_bonded_to_P = true;
			if ( distN2 < distN1 ) {
				N1_bonded_to_P = false;
			}
			id::AtomID N_atom_id;
			if ( N1_bonded_to_P ) {
				N_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("NH1") , pR_resnum );
			} else {
				N_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("NH2") , pR_resnum );
			}

			min_pose.add_constraint( new core::scoring::constraints::CoordinateConstraint( N_atom_id, fixed_atom_id, pTyr_OH.xyz(), func ) );

			TR << "Running first minimization for pArg rotamer " << model << "...." << std::endl;

			/*
			{ // scope constraint on O1P of pTyr
				// figure out which pArg oxygen is closest to O1P of pTyr
				id::AtomID O_atom_id;
				core::Real const dist1 = distance(min_pose.residue(pR_resnum).atom( "O1P" ).xyz(),pTyr_O1P.xyz() );
				core::Real const dist2 = distance(min_pose.residue(pR_resnum).atom( "O2P" ).xyz(),pTyr_O1P.xyz() );
				core::Real const dist3 = distance(min_pose.residue(pR_resnum).atom( "O3P" ).xyz(),pTyr_O1P.xyz() );
				if ( dist3 < dist1 ) {
					if ( dist2 < dist3 ) {
						O_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("O2P") , pR_resnum );
					} else {
						O_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("O3P") , pR_resnum );
					}
				} else {
					if ( dist2 < dist1 ) {
						O_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("O2P") , pR_resnum );
					} else {
						O_atom_id = id::AtomID( min_pose.residue(pR_resnum).atom_index("O1P") , pR_resnum );
					}
				}
				min_pose.add_constraint( new core::scoring::constraints::CoordinateConstraint( O_atom_id, fixed_atom_id, pTyr_O1P.xyz(), func ) );
			} // scope constraint on O1P of pTyr
			*/

			core::optimization::AtomTreeMinimizer minimizer;
			core::optimization::MinimizerOptions min_options( "dfpmin", 0.00001, true, false );
			minimizer.run( min_pose, mm_pArg, *scorefxn, min_options );

			// Throw out those where the desired P is off by > 0.75 A or the desired N is off by > 1.5 A
			if ( distance(min_pose.residue(pR_resnum).atom( "P" ).xyz(),pTyr_P_xyz) > 0.75 ) continue;

			if ( N1_bonded_to_P ) {
				if ( distance(min_pose.residue(pR_resnum).atom( "NH1" ).xyz(),pTyr_OH.xyz()) > 1.5 ) continue;
			} else {
				if ( distance(min_pose.residue(pR_resnum).atom( "NH2" ).xyz(),pTyr_OH.xyz()) > 1.5 ) continue;
			}

			// JK DEBUG
			min_pose.dump_scored_pdb( basic::options::start_file()+".int1min"+utility::to_string(model)+".pdb", *scorefxn );

			TR << "Running second minimization for pArg rotamer " << model << "...." << std::endl;
			//			min_pose.remove_constraints();
			kinematics::MoveMap mm_all;
			mm_all.set_chi( true );
			mm_all.set_bb( false );
			mm_all.set_jump( false );
			minimizer.run( min_pose, mm_all, *scorefxn, min_options );

			// JK DEBUG
			min_pose.dump_scored_pdb( basic::options::start_file()+".int2min"+utility::to_string(model)+".pdb", *scorefxn );

			TR << "Running third minimization for pArg rotamer " << model << "...." << std::endl;
			min_pose.remove_constraints();
			minimizer.run( min_pose, mm_all, *scorefxn, min_options );

			// dump scored output pdb
			TR << "Printing minimized structure " << model << std::endl;
			min_pose.dump_scored_pdb( basic::options::start_file()+".aftermin"+utility::to_string(model)+".pdb", *scorefxn );

		}

	}
	out_stream.close();
	TR << "Found " << model << " good Arg rotamers" << std::endl;

	TR << "Successfully finished building Args" << std::endl;

	return 0;
}


