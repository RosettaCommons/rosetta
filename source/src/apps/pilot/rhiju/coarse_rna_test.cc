// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief


// libRosetta headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>


//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <protocols/idealize/idealize.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/loops/Loop.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>

#include <protocols/toolbox/AllowInsert.hh>
#include <protocols/farna/MultipleDomainMover.hh>
#include <protocols/farna/RNA_ChunkLibrary.hh>
#include <protocols/farna/RNA_DataReader.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>
#include <protocols/farna/RNA_SecStructInfo.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/stepwise/enumerate/rna/StepWiseRNA_Util.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/coarse_rna/CoarseRNA_DeNovoProtocol.hh>
#include <protocols/coarse_rna/CoarseRNA_LoopCloser.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
//#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <basic/Tracer.hh>
using basic::T;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using ObjexxFCL::format::A;
using numeric::conversions::degrees;

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key

OPT_KEY( Boolean, frag_test )
OPT_KEY( Boolean, icoor )
OPT_KEY( Boolean, pdbstats )
OPT_KEY( Boolean, convert )
OPT_KEY( Boolean, jump_database )
OPT_KEY( Boolean, output_coords )
OPT_KEY( Boolean, no_sim_anneal )
OPT_KEY( Boolean, staged_constraints )
OPT_KEY( Boolean, close_loop_test )
OPT_KEY( Boolean, close_loops )
OPT_KEY( Boolean, choose_best_solution )
OPT_KEY( Boolean, force_ideal_chainbreak )
OPT_KEY( Boolean, check_pairing_dists )
OPT_KEY( Boolean, skip_base_pair_constraints )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, rb_test )
OPT_KEY( Boolean, enumerate )
OPT_KEY( Boolean, sample_angles )
OPT_KEY( Boolean, output_vdw_pose )
OPT_KEY( Boolean, little_motif )
OPT_KEY( Boolean, tar_motif )
OPT_KEY( Boolean, mismatch )
OPT_KEY( Boolean, freeze_domains )
OPT_KEY( Boolean, coarse_to_full )
OPT_KEY( Real, temperature )
OPT_KEY( Real, bin_width )
OPT_KEY( Integer, cycles )
OPT_KEY( Integer, nbulge )
OPT_KEY( String, params_file )
OPT_KEY( String,  data_file )
OPT_KEY( String,  cst_file )
OPT_KEY( IntegerVector, input_res )

using core::pose::ResMap;

/////////////////////////////////////////////////
void
coarse_frag_test(){

 	using namespace core::scoring;
 	using namespace core::chemical;
	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::io::silent;

	// read in desired sequence.
	core::sequence::Sequence fasta_sequence = *(core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]);

	// create extended coarse grained pose.
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "coarse_rna" );
	Pose pose;
	protocols::farna::make_extended_coarse_pose( pose, fasta_sequence.sequence() );

	// visualize it.
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	// (Later) --> try to set up with structure parameters.
	// read in source of fragments
	std::string infile  = option[ in ::file::s ][1];
	Pose frag_source_pose;
	import_pose::pose_from_pdb( frag_source_pose, *rsd_set, infile );

	// start doing some copy_dofs at an arbitrary position
	Size insert_res( 4 );
	for( Size n = 1; n <= frag_source_pose.total_residue()-2 ; n++ ){

		std::map< AtomID, AtomID > atom_id_map;
		std::map< Size, Size > res_map;

		for ( Size offset = 0; offset <= 2; offset ++ ) {

			res_map[ insert_res + offset ] =  n + offset;

			for ( Size j = 1; j <= pose.residue_type( insert_res + offset ).natoms(); j++ ) {
				if ( offset  == 2 && j > 1 ) continue; //only want phosphate from third residue.
				atom_id_map[  AtomID( j, insert_res + offset) ] = AtomID( j, n + offset );
			}

		}

		copy_dofs(  pose, frag_source_pose, atom_id_map );

		if ( n >= 10 && n <= 15 ) pose.dump_pdb( "S_"+string_of(  n) +".pdb" );

	}

	// score it.

}


/////////////////////////////////////////////////
void
pdb_stats( pose::Pose const & pose, utility::vector1< protocols::stepwise::enumerate::rna::PuckerState > const & pucker_states ) {

	using namespace core::id;
	using namespace core::chemical::rna;
	using namespace core::kinematics;
	using numeric::angle_radians;
	using numeric::conversions::degrees;
	using namespace protocols::stepwise::enumerate::rna;

	utility::io::ozstream out1( "dists_angles_torsions.txt" );
	utility::io::ozstream out2( "cen_xyz_diff_frames.txt" );

	Real dist_P_S( 0.0 ), dist_S_P( 0.0 ), dist_S_CEN( 0.0 ), dist_P_P( 0.0 ), dist_S_S( 0.0 );
	Real angle_S_P_S( 0.0 ), angle_P_S_P( 0.0 ), angle_P_S_CEN( 0.0 ), angle_CEN_S_P( 0.0 );
	Real torsion_S_P_S_P( 0.0 ), torsion_P_S_P_S( 0.0 );
	PuckerState prev_pucker( NORTH ), current_pucker( NORTH ), next_pucker( NORTH );

	for ( Size n = 1; n <= pose.total_residue(); n++ ) {

		if ( n > 1 ) prev_pucker = pucker_states[ n-1];
		current_pucker = pucker_states[ n ];
		if ( n < pose.total_residue() ) next_pucker = pucker_states[ n+1 ];

		dist_P_S = ( pose.xyz( NamedAtomID( " P  ", n ) ) - pose.xyz( NamedAtomID( " S  ", n ) ) ).length();
		dist_S_CEN = ( pose.xyz( NamedAtomID( " S  ", n ) ) - pose.xyz( NamedAtomID( " CEN", n ) ) ).length();

		if ( n < pose.total_residue() ) {
			dist_S_P = ( pose.xyz( NamedAtomID( " S  ", n ) ) - pose.xyz( NamedAtomID( " P  ", n+1 ) ) ).length();
			dist_P_P = ( pose.xyz( NamedAtomID( " P  ", n ) ) - pose.xyz( NamedAtomID( " P  ", n+1 ) ) ).length();
			dist_S_S = ( pose.xyz( NamedAtomID( " S  ", n ) ) - pose.xyz( NamedAtomID( " S  ", n+1 ) ) ).length();
		}

		angle_P_S_CEN = degrees(angle_radians( pose.xyz( NamedAtomID( " P  ", n) ),
																					 pose.xyz( NamedAtomID( " S  ", n) ),
																					 pose.xyz( NamedAtomID( " CEN", n) ) ) );
		if ( n < pose.total_residue() ) {
			angle_S_P_S = degrees(angle_radians( pose.xyz( NamedAtomID( " S  ", n) ),
																					 pose.xyz( NamedAtomID( " P  ", n+1) ),
																					 pose.xyz( NamedAtomID( " S  ", n+1) ) ) );
			angle_P_S_P = degrees(angle_radians( pose.xyz( NamedAtomID( " P  ", n) ),
																					 pose.xyz( NamedAtomID( " S  ", n) ),
																					 pose.xyz( NamedAtomID( " P  ", n+1) ) ) );
			angle_CEN_S_P = degrees(angle_radians( pose.xyz( NamedAtomID( " CEN", n) ),
																						 pose.xyz( NamedAtomID( " S  ", n) ),
																						 pose.xyz( NamedAtomID( " P  ", n+1) ) ) );
		}


		if ( n < pose.total_residue() ) {
			torsion_P_S_P_S = dihedral_degrees( pose.xyz( NamedAtomID( " P  ", n) ),
																					pose.xyz( NamedAtomID( " S  ", n) ),
																					pose.xyz( NamedAtomID( " P  ", n+1) ),
																					pose.xyz( NamedAtomID( " S  ", n+1) ) );
		}
		if ( n < pose.total_residue()-1 ) {
			torsion_S_P_S_P = dihedral_degrees(
																					pose.xyz( NamedAtomID( " S  ", n) ),
																					pose.xyz( NamedAtomID( " P  ", n+1) ),
																					pose.xyz( NamedAtomID( " S  ", n+1) ),
																					pose.xyz( NamedAtomID( " P  ", n+2) ) );
		}

		out1 << is_purine(  pose.residue( n ) );
		out1 << ' ' << prev_pucker <<  ' ' << current_pucker << ' ' << next_pucker;
		out1 << ' ' << dist_P_S << ' ' << dist_S_P << ' ' << dist_S_CEN << ' ' << dist_P_P << ' ' << dist_S_S;
		out1 << ' ' << angle_S_P_S << ' ' << angle_P_S_P << ' ' << angle_P_S_CEN << ' ' << angle_CEN_S_P;
		out1 << ' ' << torsion_P_S_P_S << ' ' << torsion_S_P_S_P;
		out1 << std::endl;

		if ( n > 1 && n < pose.total_residue()) {
			Stub stub1(	pose.xyz( NamedAtomID( " S  ", n) ),
									pose.xyz( NamedAtomID( " S  ", n) ),
									pose.xyz( NamedAtomID( " P  ", n) ),
									pose.xyz( NamedAtomID( " S  ", n-1) )
									);
			Vector cen_vec1 = stub1.global2local( pose.xyz( NamedAtomID( " CEN", n) ) );

			Stub stub2(	pose.xyz( NamedAtomID( " S  ", n) ),
									pose.xyz( NamedAtomID( " S  ", n) ),
									pose.xyz( NamedAtomID( " P  ", n) ),
									pose.xyz( NamedAtomID( " P  ", n+1) )
									);
			Vector cen_vec2 = stub2.global2local( pose.xyz( NamedAtomID( " CEN", n) ) );

			out2 << is_purine( pose.residue( n ) ) << ' '
					 << cen_vec1( 1 ) << ' ' <<  cen_vec1( 2 ) << ' ' << cen_vec1( 3 ) <<
				' ' << cen_vec2( 1 ) << ' ' <<  cen_vec2( 2 ) << ' ' << cen_vec2( 3 )  << std::endl;
		}

	}

	out1.close();
	out2.close();

}

/////////////////////////////////////////////////
void
vdw_stats( pose::Pose const & pose ) {

	using namespace core::conformation;
	using namespace core::id;
	using namespace core::chemical::rna;
	using namespace core::kinematics;

	utility::io::ozstream out1( "vdw.txt" );
	Real const CUTOFF( 12.0 );

	Size const NATOMS( 3 ); //P, S, CEN
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		Residue const & rsd1 = pose.residue( i );

		for ( Size j = i+2; j <= pose.total_residue(); j++ ) {
			Residue const & rsd2 = pose.residue( j );

			utility::vector1< Real > dists;
			Real min_dist( CUTOFF );
			for ( Size k = 1; k <= NATOMS; k++ ) {
				for ( Size m = 1; m <= NATOMS; m++ ) {

					Real const dist = ( rsd1.xyz( k ) - rsd2.xyz( m ) ).length();
					dists.push_back( dist );
					if ( dist < min_dist ) min_dist = dist;

				}
			}

			if ( min_dist < CUTOFF ){
				for ( Size k = 1; k <= dists.size(); k++ ) {
					out1 << ' ' << dists[ k ];
				}
				out1 << std::endl;
			}

		}
	}

	out1.close();

}

/////////////////////////////////////////////////
void
icoor_test(){
 	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose, coarse_pose;

	////////////////////////////////
	// Figure out icoor
	std::string full_sequence = "acguacguacgu";
	make_pose_from_sequence( pose, full_sequence, *rsd_set );
	protocols::farna::make_coarse_pose( pose, coarse_pose );
	coarse_pose.dump_pdb( "extended.pdb" );

	protocols::farna::make_extended_coarse_pose( coarse_pose, pose.sequence() );
	coarse_pose.dump_pdb( "coarse_extended.pdb" );
}

/////////////////////////////////////////////////
void
convert_to_coarse_test(){
 	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace protocols::farna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose, coarse_pose;

	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );
	figure_out_reasonable_rna_fold_tree( pose );

	protocols::farna::make_coarse_pose( pose, coarse_pose );
	std::cout << "--------------------" << std::endl;
	std::cout << "Check out coarse.pdb" << std::endl;
	std::cout << "--------------------" << std::endl;
	coarse_pose.dump_pdb( "coarse.pdb" );


	SilentFileData silent_file_data;

	BinaryRNASilentStruct s1( pose, "S_0" );
	silent_file_data.write_silent_struct( s1, "fullatom.out", false );

	BinaryRNASilentStruct s2( coarse_pose, "S_0" );
	silent_file_data.write_silent_struct( s2, "coarse.out", false );

}


/////////////////////////////////////////////////
void
output_minipose_coords_test(){
 	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	protocols::farna::figure_out_secstruct( pose );
	std::string secstruct( protocols::farna::get_rna_secstruct( pose ) );

	utility::io::ozstream out1( "coarse_coords.txt" );

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		out1 << pose.sequence()[i-1] << " " << secstruct[ i-1 ];
		for ( Size j = 1; j <= pose.residue(i).natoms(); j++ ) {
			Vector coord = pose.residue(i).xyz(j);
			out1 << F(12,6,coord(1)) << F(12,6,coord(2)) << F(12,6,coord(3));
		}
		out1 << std::endl;
	}


	std::cout << "---------------------------------------" << std::endl;
	std::cout << "Output coordinates to coarse_coords.txt" << std::endl;
	std::cout << "---------------------------------------" << std::endl;



}

/////////////////////////////////////////////////
void
pdbstats_test(){

 	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace protocols::stepwise::enumerate::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( RNA );

	////////////////////////////////
	Pose pose, coarse_pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	utility::vector1< PuckerState > pucker_states;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) pucker_states.push_back( Get_residue_pucker_state( pose, n ) );

	protocols::farna::make_coarse_pose( pose, coarse_pose );
	coarse_pose.dump_pdb( "coarse.pdb" );

	// How about some pdb_stats?
	pdb_stats( coarse_pose, pucker_states );
	vdw_stats( coarse_pose );

}

//////////////////////////////////////////////////////////////////////////////////////
// JUMP extractor.
void
create_bp_jump_database_test( ){

	using namespace chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	std::string infile  = option[ in::file::s ][1];
	std::string outfile  = option[ out::file::o ];

	pose::Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	// Fill base pairing information... these are
	// all functions used in scoring... see RNA_BaseBaseEnergy.cc
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	utility::io::ozstream dataout( outfile );

	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		int const i = base_pair.res1;
		int const k = base_pair.edge1;

		int const j = base_pair.res2;
		int const m = base_pair.edge2;

		char const orientation = ( base_pair.orientation == 1) ? 'A' : 'P';

		char const edge_i = get_edge_from_num( k );
		char const edge_j = get_edge_from_num( m );

		//Figure out jump.
		conformation::Residue const & rsd_i( pose.residue( i ) );
		conformation::Residue const & rsd_j( pose.residue( j ) );
		kinematics::Stub const stub_i( rsd_i.xyz( " Y  " ),
																	 rsd_i.xyz( " X  " ),
																	 rsd_i.xyz( " CEN" ) );
		kinematics::Stub const stub_j( rsd_j.xyz( " Y  " ),
																	 rsd_j.xyz( " X  " ),
																	 rsd_j.xyz( " CEN" ) );

		dataout << "PAIR " <<
			I(5, i) << ' ' << edge_i << ' ' <<
			I(5, j) << ' ' << edge_j << "   " <<
			orientation << "   " <<
			pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
			" Y  "  <<  " " <<
			" Y  "  <<  " " <<
		  kinematics::Jump( stub_i, stub_j) <<
			std::endl;

	}

	dataout.close();

	std::cout << "***********************************************************" << std::endl;
	std::cout << "Put jumps from PDB file " <<  infile << " into " << outfile << std::endl;
	std::cout << "***********************************************************" << std::endl;

}


//////////////////////////////////////////////////////////////////////////////////////
void
general_initialize( 	pose::Pose & pose,
											pose::PoseOP & native_pose,
											protocols::farna::RNA_StructureParametersOP & rna_structure_parameters_,
											protocols::coarse_rna::CoarseRNA_LoopCloserOP & rna_loop_closer_,
											protocols::farna::RNA_ChunkLibraryOP & rna_chunk_library_,
											protocols::toolbox::AllowInsertOP &  allow_insert_
){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace core::pose;

	///////////////////////////////////////////////////////////////////////////
	// pdb setup
	///////////////////////////////////////////////////////////////////////////
	// read in desired sequence.
	core::sequence::Sequence fasta_sequence = *(core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]);

	protocols::farna::make_extended_coarse_pose( pose, fasta_sequence.sequence() );
	pose.dump_pdb( "extended.pdb" );

	// native?
	if ( option[ in::file::native ].user() ) {
		ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
		std::string native_pdb_file  = option[ in::file::native ];
		native_pose = new pose::Pose;
		import_pose::pose_from_pdb( *native_pose, *rsd_set_coarse, native_pdb_file );
	}

	///////////////////////////////////////////////////////////////////////////
	// use parameters, copy dofs for pose setup.
	///////////////////////////////////////////////////////////////////////////

	std::string const rna_params_file_( option[ params_file ]() );
	std::string const jump_library_file_( basic::database::full_name("sampling/rna/1jj2_coarse_jumps.dat" ) );
	utility::vector1< std::string > chunk_silent_files_( option[ in::file::silent ]() );

	rna_structure_parameters_ = new RNA_StructureParameters;
	rna_structure_parameters_->initialize( pose, rna_params_file_, jump_library_file_, true /*ignore_secstruct*/ );

	utility::vector1< Size > input_res_( option[ input_res ]() );
	if( input_res_.size() > 0 ) {
		rna_chunk_library_ = new protocols::farna::RNA_ChunkLibrary( chunk_silent_files_, pose, input_res_ );
	} else {
		rna_chunk_library_ = new protocols::farna::RNA_ChunkLibrary( chunk_silent_files_, pose, rna_structure_parameters_->connections() );
	}

	rna_structure_parameters_->set_allow_insert( rna_chunk_library_->allow_insert() );
	rna_structure_parameters_->setup_fold_tree_and_jumps_and_variants( pose );

	pose.dump_pdb( "after_fold_tree.pdb" );

	rna_chunk_library_->initialize_random_chunks( pose, true /*dump_pdb*/ );

	std::cout << "ALLOW INSERT " << std::endl;
	allow_insert_ = rna_structure_parameters_->allow_insert();
	allow_insert_->show();
	rna_loop_closer_ = new CoarseRNA_LoopCloser;
	rna_loop_closer_->set_allow_insert( allow_insert_ );

	//	if ( choose_best_solution_ ) rna_loop_closer_->choose_best_solution_based_on_score_function( denovo_scorefxn_ );

	// read in helix PDB
	pose.dump_pdb( "start.pdb" );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

}



//////////////////////////////////////////////////////////////////////////////////////
// JUMP extractor.
void
coarse_rb_test(){

	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace protocols::moves;
	using namespace protocols::toolbox;
	using namespace core::pose;

	// create extended coarse grained pose.
	pose::Pose pose;
	pose::PoseOP native_pose;
	RNA_StructureParametersOP rna_structure_parameters_;
	CoarseRNA_LoopCloserOP rna_loop_closer_;
	RNA_ChunkLibraryOP rna_chunk_library_;
	AllowInsertOP allow_insert_;

	general_initialize( pose, native_pose, rna_structure_parameters_, rna_loop_closer_, rna_chunk_library_, allow_insert_ );

	MultipleDomainMoverOP multiple_domain_mover = new MultipleDomainMover( pose, rna_loop_closer_ );
	Size const num_domains = multiple_domain_mover->num_domains();

	multiple_domain_mover->randomize_pose_rigid_bodies( pose );

	// do rigid body moves to 1, rigid body moves to 2.
	// close loops?
	Size const cycles( 10000 );
	Size const iter( 1000);
	for ( Size i = 1; i <= iter; i++ ) {

		for ( Size n = 1; n <= num_domains; n++ ) {

			std::cout << "DOING DOMAIN " << n << std::endl;

			for ( Size j = 1; j <= cycles; j++ ) {
				Size const jumpno = multiple_domain_mover->apply_at_domain( pose, n );
				rna_loop_closer_->apply_after_jump_change( pose, jumpno );
			}

		}

	}

	pose.dump_pdb( "final.pdb" );


}

//////////////////////////////////////////////////////////////////////////
void
output_angles( core::pose::Pose const & pose ){
	using namespace numeric;
	using namespace numeric::conversions;
	using namespace core::id;
	using namespace core;

	core::Real theta_angle;
	std::cout << "THETA ";

	for ( core::Size i = 1; i < pose.total_residue(); i++ ) {
		theta_angle = degrees( angle_radians( pose.xyz( NamedAtomID( " P  ", i) ),
																		pose.xyz( NamedAtomID( " S  ", i) ),
																		pose.xyz( NamedAtomID( " P  ", i+1) ) ) );
		std::cout << ' ' << theta_angle;

		theta_angle = degrees( angle_radians( pose.xyz( NamedAtomID( " S  ", i) ),
																		pose.xyz( NamedAtomID( " P  ", i+1) ),
																		pose.xyz( NamedAtomID( " S  ", i+1) ) ) );
		std::cout << ' ' << theta_angle;
	}

	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////
void
get_ideal_angle_sets(  	 utility::vector1< utility::vector1< Real > > & angle_sets,
												 utility::vector1< utility::vector1< core::id::AtomID > > const & angle_ids,
												 utility::vector1< Real >  angle_set ){

	using numeric::conversions::radians;

	if ( angle_set.size() == angle_ids.size() ) {
		angle_sets.push_back( angle_set );
	} else {

		angle_set.push_back( radians( 87.0 ) );
		get_ideal_angle_sets( angle_sets, angle_ids, angle_set );

		angle_set[ angle_set.size() ] = radians( 110.0 );
		get_ideal_angle_sets( angle_sets, angle_ids, angle_set );

	}

}

//////////////////////////////////////////////////////////////////////////////////////
void
get_angle_sets(  pose::Pose const & pose,
								 protocols::toolbox::AllowInsertOP allow_insert,
								 utility::vector1< utility::vector1< Real > > & angle_sets,
								 utility::vector1< utility::vector1< core::id::AtomID > > & angle_ids,
								 bool const sample_angles_ ){

	using namespace core::id;

	for ( Size i = 1; i < pose.total_residue(); i++ ){

		utility::vector1< AtomID > ids;
		// P-S-P
		if ( allow_insert->get( AtomID(2,i) /*sugar*/ )  ) {
			ids.clear();
			if ( pose.fold_tree().is_cutpoint( i ) && pose.residue_type( i ).has_variant_type( "CUTPOINT_LOWER" ) ) {
				ids.push_back( AtomID( 1, i ) );
				ids.push_back( AtomID( 2, i ) );
				ids.push_back( AtomID( 3, i ) ); //OVL1
			} else {
				ids.push_back( AtomID( 1, i ) );
				ids.push_back( AtomID( 2, i ) );
				ids.push_back( AtomID( 1, i+1 ) );
			}
			angle_ids.push_back( ids );
		}

		// S-P-S
		if ( allow_insert->get( AtomID(1,i+1) /*phosphate*/ )  ) {
			ids.clear();
			if ( pose.fold_tree().is_cutpoint( i ) && pose.residue_type( i ).has_variant_type( "CUTPOINT_LOWER" ) ) {
				ids.push_back( AtomID( 2, i ) );
				ids.push_back( AtomID( 3, i ) );
				ids.push_back( AtomID( 4, i ) );
			} else {
				ids.push_back( AtomID( 2, i ) );
				ids.push_back( AtomID( 1, i+1 ) );
				ids.push_back( AtomID( 2, i+1 ) );
			}
			angle_ids.push_back( ids );
		}

	}

	std::cout << "NUM ANGLE IDS? " << angle_ids.size() << std::endl;

	utility::vector1< Real > angle_set;
	if ( sample_angles_ )  {
		get_ideal_angle_sets( angle_sets, angle_ids, angle_set );
	} else {
		for ( Size i = 1; i <= angle_ids.size(); i++ ) {
			//std::cout << angle_ids[i].size() << std::endl;

			Real angle_val = pose.conformation().bond_angle( angle_ids[i][1], angle_ids[i][2], angle_ids[i][3] );
			std::cout << angle_val << std::endl;
			angle_set.push_back( angle_val );
		}
		angle_sets.push_back( angle_set );
	}
	std::cout << "NUM ANGLE SETS " << angle_sets.size() << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
void
apply_angle_set( pose::Pose & pose,
								 utility::vector1< utility::vector1< core::id::AtomID > > const & angle_ids,
								 utility::vector1< Real > const & angle_set ){
	using namespace core::id;
	std::cout << "ANGLE_SET_SIZE " << angle_set.size() << " " << angle_ids.size() << std::endl;

	for ( Size i = 1; i <= angle_set.size(); i++ ) {
		std::cout << "ANGLE: " << i << std::endl;

		pose.conformation().set_bond_angle( angle_ids[i][1], angle_ids[i][2], angle_ids[i][3], angle_set[ i ] );

		// FOLLOWING DOES NOT WORK; WOULD NEED TO USE SPHERICAL COORDS CRAP.
		//		if ( angle_ids[i][3].atomno() == 4 ) {  //special case --> OVL2. Also change next guy's bond angle
		//			Size const n = angle_ids[i][3].rsd();
		//			std::cout << "SPECIAL CASE " << n << pose.residue_type( n+1).has_variant_type( chemical::CUTPOINT_UPPER ) << std::endl;
		//			pose.conformation().set_bond_angle( AtomID( n+1, 2 ) /*S=OVU1*/, AtomID( n+1, 1 ) /*P*/, AtomID( n+1, 3 ) /*S*/ , angle_set[ i ] );
		//		}

	}
}

//////////////////////////////////////////////////////////////////////////////////////
void
output_torsions( pose::Pose const & coarse_pose, Size const & i, utility::io::ozstream & out1 ){

	using namespace core::id;
	using numeric::conversions::degrees;

	Real torsion1 =  coarse_pose.torsion( TorsionID( i, BB, 2 ) );
	if ( torsion1 < 0.0 ) torsion1 += 360.0;

	Real torsion2 =  coarse_pose.torsion( TorsionID( i+1, BB, 1 ) );
	if ( torsion2 < 0.0 ) torsion2 += 360.0;

	out1 << torsion1 << ' ' << torsion2 << ' ';

	out1 << ' ' << degrees( coarse_pose.conformation().bond_angle( AtomID(1,i), AtomID(2,i), AtomID(1,i+1) ) );
	out1 << ' ' << degrees( coarse_pose.conformation().bond_angle( AtomID(2,i), AtomID(1,i+1), AtomID(2,i+1) ) );

	out1 << ' ' << i << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////
void
pdb_little_motif_test(){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace core::pose;

	// native?
	pose::Pose pose, coarse_pose;

	ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( RNA );
	std::string infile  = option[ in::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set_coarse, infile );

	protocols::farna::make_coarse_pose( pose, coarse_pose );

	std::map< Size, Size > partner;
	figure_out_base_pair_partner( pose, partner );


	utility::io::ozstream out( "all_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-1; i++ ) {
		output_torsions( coarse_pose, i, out );
	}

	// Two-helix junction, 1x1.
	utility::io::ozstream out1( "two_way_1_1_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-2; i++ ) {
		if ( partner.find(i  ) == partner.end() ) continue;
		std::cout << i << " " <<  partner[ i ] << std::endl;
		if ( partner.find(i+2) == partner.end() ) continue;
		if ( partner[ i ] != partner[ i+2 ] +2 ) continue;

		output_torsions( coarse_pose, i, out1);

	}
	out1.close();


	// Two-helix junction, 1x1.
	utility::io::ozstream out2( "two_way_0_2_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-2; i++ ) {

		if ( partner.find(i  ) == partner.end() ) continue;
		if ( partner.find(i+1) == partner.end() ) continue;
		if ( partner[ i ] != partner[ i+1 ] +3 ) continue;

		output_torsions( coarse_pose, i, out2);

	}
	out2.close();

	// Triloop.
	utility::io::ozstream out3( "triloop_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-4; i++ ) {
		if ( partner.find(i  ) == partner.end() ) continue;
		if ( partner[ i ] != i+4 ) continue;

		output_torsions( coarse_pose, i, out3);

	}
	out3.close();

	// Three-way junction
	utility::io::ozstream out4( "three_way_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-1; i++ ) {
		if ( partner.find( i    ) == partner.end() ) continue;
		if ( partner.find( i+1  ) == partner.end() ) continue;
		if ( partner.find( partner[i+1] + 2  ) == partner.end() ) continue;
		if ( partner.find( partner[ partner[i+1] + 2] + 1 ) == partner.end() ) continue;
		if ( partner[partner[partner[i+1] + 2] + 1 ] != i ) continue;

		output_torsions( coarse_pose, i, out4);
	}
	out4.close();

	// Four-way junction
	utility::io::ozstream out5( "four_way_STATS.txt" );
	for ( Size i = 1; i <= pose.total_residue()-2; i++ ) {
		if ( partner.find( i    ) == partner.end() ) continue;
		if ( partner.find( i+1  ) == partner.end() ) continue;
		if ( partner.find( partner[i+1] +1 ) == partner.end() ) continue;
		if ( partner.find( partner[ partner[i+1] +1 ] + 1) == partner.end() ) continue;
		if ( partner.find( partner[ partner[ partner[i+1] +1 ] + 1 ] + 1 ) == partner.end() ) continue;
		if ( partner[ partner[ partner[ partner[i+1] +1 ] + 1 ] + 1 ] != i ) continue;

		output_torsions( coarse_pose, i, out5);
	}
	out5.close();

}


//////////////////////////////////////////////////////////////////////////////////////
void
get_base_pair_coordinate_system( pose::Pose const & pose, Size const res1, Size const res2, Vector & centroid, Matrix & M ){

	using namespace core::id;
	using namespace core::chemical;

	Vector x = pose.xyz( NamedAtomID(" C1'", res2 ) ) -  pose.xyz( NamedAtomID(" C1'", res1 )  );
	x.normalize();

	centroid = 0.5 * ( pose.xyz( NamedAtomID(" C1'", res2 ) ) +  pose.xyz( NamedAtomID(" C1'", res1 )  ) );

  Size res_type = pose.residue( res1 ).aa();

	// Make a perpendicular axis pointing from centroid towards
	// Hoogstein edge (e.g., major groove in a double helix).
	std::string H_atom;
	if ( res_type == na_rad ) H_atom = "N7";
	if ( res_type == na_rcy ) H_atom = "C5";
	if ( res_type == na_rgu ) H_atom = "N7";
	if ( res_type == na_ura ) H_atom = "C5";

	Vector y = pose.residue( res1 ).xyz( H_atom ) - centroid;
	Vector z = cross( x, y );
	z.normalize();

	y = cross( z, x );

	for ( Size k = 1; k <= 3; k++ ){
		M(k,1) =x(k);
		M(k,2) =y(k);
		M(k,3) =z(k);
	}


}

//////////////////////////////////////////////////////////////////////////////////////
void
reorient_to_base_pair_coordinate_system( pose::Pose & pose, Size const res1, Size const res2 ){

	using namespace core::id;

	Vector centroid;
	Matrix M;
	get_base_pair_coordinate_system( pose, res1, res2, centroid, M );

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ){
			pose.set_xyz(  AtomID( j,i ),
										 M.transposed() * (pose.xyz( AtomID( j, i ) ) - centroid) );
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////
Matrix
cycle( Matrix const & M0 ){
	Matrix M;
	for ( Size k = 1; k <= 3; k++ ){
		M(k,1) = M0(k,3);
		M(k,2) = M0(k,2);
		M(k,3) = M0(k,1);
	}
	return M;
}

//////////////////////////////////////////////////////////////////////////////////////
void
determine_delx_dely_delz_alpha_beta_gamma( pose::Pose const & pose,
																					 Size const bp1_partner1,
																					 Size const bp1_partner2,
																					 Size const bp2_partner1,
																					 Size const bp2_partner2,
																					 utility::io::ozstream & out,
																					 bool const cycle_xyz_to_yzx = false ){

	Vector centroid1, centroid2;
	Matrix M1, M2;
	get_base_pair_coordinate_system( pose, bp1_partner1, bp1_partner2, centroid1, M1 );
	get_base_pair_coordinate_system( pose, bp2_partner1, bp2_partner2, centroid2, M2 );

	if (cycle_xyz_to_yzx){
		M1 = cycle( M1 );
		M2 = cycle( M2 );
	}

	Real const	delx = centroid2(1) - centroid1(1);
	Real const	dely = centroid2(2) - centroid1(2);
	Real const	delz = centroid2(3) - centroid1(3);

	Real alpha, beta, gamma;
	protocols::stepwise::get_euler_angles( alpha, beta, gamma, M1, M2, false /*verbose*/ );

	out << alpha << ' ' << beta << ' ' << gamma << ' ' << delx << ' ' << dely << ' ' << delz << std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////
void
tar_motif_test(){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace core::pose;

	pose::Pose pose;

	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( RNA );
	utility::vector1< std::string > const & infiles  = option[ in::file::s ]();

	Size const NBULGE = option[ nbulge ]();


	// Four-way junction
	std::string const outfile = option[ out::file::o ]();
	std::string const outfile1 = "xyz_" + outfile;
	std::string const outfile2 = "yzx_" + outfile;
	utility::io::ozstream out1( outfile1 );
	utility::io::ozstream out2( outfile2 );

	for( Size n = 1; n <= infiles.size(); n++ ) {

		std::string infile  = infiles[n];
		std::cout << " about to read from: " << infile << std::endl;

		import_pose::pose_from_pdb( pose, *rsd_set, infile );

		std::map< Size, Size > partner;

		if ( pose.total_residue() < 3 ) continue;

		figure_out_base_pair_partner( pose, partner );

		// look for three residue bulges.
		for ( Size i = 1; i <= pose.total_residue()-1; i++ ) {

			utility::vector1< Size > motif_res;

			if ( partner.find( i    ) == partner.end() ) continue;
			if ( partner.find( i+1  ) == partner.end() ) continue;

			motif_res.push_back( i );
			motif_res.push_back( i+1 );

			if ( partner.find( partner[i+1] + NBULGE + 1  ) == partner.end() ) continue;
			if ( partner[ partner[i+1] + NBULGE + 1 ] != i ) continue;

			for ( Size k = 0; k <= NBULGE+1; k++ ){
				motif_res.push_back( partner[i+1] + k );
			}

			Pose mini_pose;
			protocols::stepwise::pdbslice( mini_pose, pose, motif_res );

			reorient_to_base_pair_coordinate_system( mini_pose, 1, motif_res.size() );

			out1 << infile << ' ' << i << ' ' << partner[ i+1 ] << ' ' ;
			determine_delx_dely_delz_alpha_beta_gamma( mini_pose, 1, motif_res.size(), 2, 3, out1 );

			out2 << infile << ' ' << i << ' ' << partner[ i+1 ] << ' ' ;
			determine_delx_dely_delz_alpha_beta_gamma( mini_pose, 1, motif_res.size(), 2, 3, out2, true /* xyz --> yzx */ );

			Size pos( infile.find( "_RNA.pdb" ) );

			std::string const filename = infile.substr( pos-4 , 4 ) + "_" +
				ObjexxFCL::lead_zero_string_of( i,4 ) + "_" +
				ObjexxFCL::lead_zero_string_of( partner[i+1], 4 ) + ".pdb";

			mini_pose.dump_pdb( filename );

		}


	}

	out1.close();
	out2.close();

}

///////////////////////////////////////////////////////////////////////////////
void
mismatch_test(){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace core::pose;

	pose::Pose pose;

	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( RNA );
	utility::vector1< std::string > const & infiles  = option[ in::file::s ]();

	//	std::string const outfile = option[ out::file::o ]();
	//	utility::io::ozstream out( outfile );

	for( Size n = 1; n <= infiles.size(); n++ ) {

		std::string infile  = infiles[n];
		std::cout << " about to read from: " << infile << std::endl;

		import_pose::pose_from_pdb( pose, *rsd_set, infile );

		std::map< Size, Size > partner;

		if ( pose.total_residue() < 6 ) continue;

		figure_out_base_pair_partner( pose, partner, false /*strict*/ );

		// look for mismatches
		for ( Size i = 1; i <= pose.total_residue()-1; i++ ) {

			utility::vector1< Size > motif_res;

			if ( partner.find( i    ) == partner.end() ) continue;
			if ( partner.find( i+1  ) != partner.end() ) continue; //force a mismatch in between
			if ( partner.find( i+2  ) == partner.end() ) continue;

			if ( partner[i+2] == i ) continue; //weird case.
			if ( is_rna_chainbreak( pose, i ) ) continue;
			if ( is_rna_chainbreak( pose, i+1 ) ) continue;

			if ( !possibly_canonical_strict( pose.aa(i),   pose.aa(partner[i]  ) ) ) continue;
			if ( !possibly_canonical_strict( pose.aa(i+2), pose.aa(partner[i+2]) ) ) continue;

			for ( Size k = 0; k <= 2; k++ )	 motif_res.push_back( i + k );

			if ( partner.find( partner[i+2] + 2  ) == partner.end() ) continue;
			if ( partner[ partner[i+2] + 2 ] != i ) continue;

			if ( is_rna_chainbreak( pose, partner[i+2] ) ) continue;
			if ( is_rna_chainbreak( pose, partner[i+1]+1 ) ) continue;

			for ( Size k = 0; k <= 2; k++ )	 motif_res.push_back( partner[i+2] + k );

			Pose mini_pose;
			protocols::stepwise::pdbslice( mini_pose, pose, motif_res );

			reorient_to_base_pair_coordinate_system( mini_pose, 1, motif_res.size() );

			Size pos( infile.find( "_RNA.pdb" ) );

			std::string const filename =
				mini_pose.sequence().substr(1,1) +
				mini_pose.sequence().substr(4,1) + "_" +
				infile.substr( pos-4 , 4 ) + "_" +
				ObjexxFCL::lead_zero_string_of( i,4 ) + "_" +
				ObjexxFCL::lead_zero_string_of( partner[i+2], 4 ) + ".pdb";
			std::cout << "Outputting mismatch pdb to: " << filename << std::endl;

			mini_pose.dump_pdb( filename );

		}


	}

}

//////////////////////////////////////////////////////////////////////////////////////
void
enumerate_map_test(){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::coarse_rna;
	using namespace protocols::farna;
	using namespace protocols::toolbox;
	using namespace core::pose;

	// create extended coarse grained pose.
	pose::Pose pose;
	pose::PoseOP native_pose;
	RNA_StructureParametersOP rna_structure_parameters_;
	CoarseRNA_LoopCloserOP rna_loop_closer_;
	RNA_ChunkLibraryOP rna_chunk_library_;
	AllowInsertOP allow_insert_;
	general_initialize( pose, native_pose, rna_structure_parameters_, rna_loop_closer_, rna_chunk_library_, allow_insert_ );

	output_angles( pose );

	// Choose first allow_insert phosphate to sample...
	Size sample_res( 0 );
	utility::vector1< TorsionID > moving_torsions;
	utility::vector1< bool > is_moving_res;

	for ( Size i = 2; i <= pose.total_residue(); i++ ){

		is_moving_res.push_back( false );

		//		bool check_domain_boundary = ( allow_insert_->get_domain( AtomID(1,i) /*phosphate*/ ) != allow_insert_->get_domain( AtomID(2,i-1) /*sugar*/) );
		//		std::cout << "CHECK " << i << ' ' << allow_insert_->get( AtomID(i,1) ) << ' ' << check_domain_boundary << std::endl;

		if// (allow_insert_->get( AtomID(i,1) ) && check_domain_boundary ){
			(allow_insert_->get( TorsionID( i-1, BB, 2 ), pose.conformation() ) &&
			 allow_insert_->get( TorsionID( i, BB, 1 ), pose.conformation() ) ){
			if ( sample_res == 0 ) sample_res = i;
			std::cout << "FOUND SAMPLE RES: " << i << std::endl;
			moving_torsions.push_back( TorsionID( i - 1, BB, 2 ) );
			moving_torsions.push_back( TorsionID( i    , BB, 1 ) );
			is_moving_res[ i ] = true;
		}

	}
	if ( sample_res > pose.total_residue() ) utility_exit_with_message( "Could not find sample_res" );
	TorsionID const & torsion1( moving_torsions[ 1 ] );
	TorsionID const & torsion2( moving_torsions[ 2 ] );

	utility::vector1< utility::vector1< AtomID > > angle_ids;
	utility::vector1< utility::vector1< Real > > angle_sets;
	bool sample_angles_( option[ sample_angles ]() );
	get_angle_sets( pose, allow_insert_, angle_sets, angle_ids, sample_angles_ );

	// virtualize base for internal residues.
	for ( Size i = 1; i < pose.total_residue(); i++ ){
		if ( is_moving_res[ i ] && is_moving_res[ i+1 ] ){
			add_variant_type_to_pose_residue( pose,  "VIRTUAL_BASE", i );
		}
	}

	if ( moving_torsions.size() != 8 ) {
		std::cout << "NUM MOVING TORSIONS " <<  moving_torsions.size() << std::endl;
		utility_exit_with_message( "Need a motif with 8 moving torsions" );
	}

	// Sample torsion angles.
	Real const bin_size( option[ bin_width]()  );
	Size const nbins(   static_cast<Size>( 360.0/ bin_size)  );
	Size count( 0 );

	//sterics.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( vdw, 1.0 );

	SilentFileData silent_file_data;
	std::string const silent_file_out_IN = option[ out::file::silent]();
	std::string silent_file_out( silent_file_out_IN );
	Real const VDW_CUTOFF( 0.1 );

	std::cout << "WILL ITERATE OVER ANGLE SETS: " << angle_sets.size() << std::endl;

	for ( Size n = 1; n <= angle_sets.size(); n++ ) {

		if ( sample_angles_ ) {
			silent_file_out = silent_file_out_IN;
			Size pos( silent_file_out.find( ".out" ) );
			silent_file_out.replace( pos, 4, "."+string_of(n)+".out" );
			apply_angle_set( pose, angle_ids, angle_sets[ n ] );
		}

		std::string const vdw_silent_file_out = "vdw_"+silent_file_out;


		for ( Size i = 1; i <= nbins; i++ ){

			Real const tau1( i * bin_size );
			pose.set_torsion( torsion1, tau1 );
			std::cout << "Doing " << I( 5,i) << " out of " <<  I( 5, nbins ) << std::endl;
			//std::cout <<  torsion1 << ' ' << tau1 << std::endl;

			for ( Size j = 1; j <= nbins; j++ ){

				Real const tau2( j * bin_size );
				pose.set_torsion( torsion2, tau2 );

				rna_loop_closer_->apply( pose, sample_res-1 );
				Size const nsol = rna_loop_closer_->nsol();

				utility::vector1< pose::PoseOP > pose_list;
				rna_loop_closer_->get_all_solutions( pose, pose_list );

				(*scorefxn)( pose ); //just to see something in the viewer.

				//std::cout << "POSE_LIST_SIZE: " << pose_list.size() << std::endl;

				for ( Size n = 1; n <= pose_list.size(); n++ ){
					Real const score = (*scorefxn)( *pose_list[n] );
					EnergyMap const & energy_map = pose_list[n]->energies().total_energies();
					Real const & vdw_score = energy_map[ vdw ];

					count++;
					std::string const tag = "S_" + lead_zero_string_of( count, 5);
					BinaryRNASilentStruct s( *pose_list[n], tag );
					s.add_energy( "nsol", nsol );


					//				output_angles( pose );

					for ( Size k = 1; k <= moving_torsions.size(); k++ ) {
						s.add_energy( "tau"+string_of( k ),  pose.torsion( moving_torsions[ k ]  ) );
					}

					if ( native_pose ){
						s.add_energy( "all_rms", all_atom_rmsd( pose, *native_pose ) );
					}

					silent_file_data.write_silent_struct( s, silent_file_out, true /*false*/ );

					if ( option[ output_vdw_pose]() && vdw_score < VDW_CUTOFF ){
						silent_file_data.write_silent_struct( s, vdw_silent_file_out, false );
					}

				}


			}
		}
	}


	// Later put in score stuff.
	//	denovo_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( lores_scorefxn_ );

}

//////////////////////////////////////////////////////////////////////////////////////
// JUMP extractor.
void
coarse_close_loop_test( ){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace protocols::coarse_rna;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	// initialize helix pose
	Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set_coarse, "ggaauucc_coarse_RNA.pdb" );

	FoldTree f( pose.total_residue() );

	f.new_jump( 4, 5, 4 );
	//	f.new_jump( 3, 6, 4 );

	Size const cutpos = 6;
	f.new_jump( 1, 8, cutpos );

	f.set_jump_atoms( 1, " Y  ", " Y  " );
	f.set_jump_atoms( 2, " Y  ", " Y  " );
	pose.fold_tree( f );

	pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpos   );
	pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpos+1 );

	pose.dump_pdb( "start.pdb" );


	// perturb pose & output.

	Size const perturb_res = 2;
	//Size const perturb_res = 5;

	pose.set_torsion( TorsionID( perturb_res,   BB, 2 ), 90.0 );
	pose.set_torsion( TorsionID( perturb_res+1, BB, 1 ), 90.0 );

	pose.dump_pdb( "perturb.pdb" );

	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( vdw, 1.0 );

	CoarseRNA_LoopCloser loop_closer;
	loop_closer.choose_best_solution_based_on_score_function( scorefxn );
	loop_closer.apply( pose, perturb_res );

	pose.dump_pdb( "closed.pdb" );

	utility::vector1< PoseOP > pose_list;
	loop_closer.get_all_solutions( pose, pose_list );
}


//////////////////////////////////////////////////////////////////////////////////////
// Totally hacky -- later would be pretty easy to work out the math
//  to convert coarse grained (two dummy atoms) to full atom, since
//  both representations have six DOFs.
void
coarse_to_full_test( ){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace protocols::coarse_rna;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	ResidueTypeSetCAP rsd_set_full = ChemicalManager::get_instance()->residue_type_set( RNA );

	// initialize helix pose
	Pose pose_coarse;
	std::string infile  = option[ in::file::s ][1];
	import_pose::pose_from_pdb( pose_coarse, *rsd_set_coarse, infile );

	Pose pose_scratch, pose_scratch_coarsened;
	make_pose_from_sequence( pose_scratch, pose_coarse.sequence(), *rsd_set_full );
	pose_scratch.dump_pdb( "scratch_full.pdb" );

	protocols::farna::make_coarse_pose( pose_scratch, pose_scratch_coarsened );
	pose_scratch.dump_pdb( "scratch_coarse.pdb" );

	Pose pose = pose_scratch;

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		// A couple coordinate systems to allow easy superposition.
		Vector const origin1 =  pose_coarse.xyz( NamedAtomID( " P  ", i ) );
		Vector const x1 =  pose_coarse.xyz( NamedAtomID( " P  ", i ) );
		Vector const y1 =  pose_coarse.xyz( NamedAtomID( " S  ", i ) );
		Vector const z1 =  pose_coarse.xyz( NamedAtomID( " CEN", i ) );
		Stub stub1( origin1, x1, y1, z1 );

		// A couple coordinate systems to allow easy superposition.
		Vector const origin2 =  pose_scratch_coarsened.xyz( NamedAtomID( " P  ", i ) );
		Vector const x2 =  pose_scratch_coarsened.xyz( NamedAtomID( " P  ", i ) );
		Vector const y2 =  pose_scratch_coarsened.xyz( NamedAtomID( " S  ", i ) );
		Vector const z2 =  pose_scratch_coarsened.xyz( NamedAtomID( " CEN", i ) );
		Stub stub2( origin2, x2, y2, z2 );

		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ){
			pose.set_xyz( AtomID(j,i), stub1.local2global( stub2.global2local( pose_scratch.xyz( AtomID(j,i) ) ) ) );
		}

	}

	std::string outfile = "full.pdb";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	std::cout<< "dumping to " << outfile << std::endl;
	pose.dump_pdb( outfile );

}



///////////////////////////////////////////////////////////////
void
coarse_rna_denovo_test(){

	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace protocols::coarse_rna;

	ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	// read in desired sequence.
	core::sequence::Sequence fasta_sequence = *(core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]);

	// create extended coarse grained pose.
	pose::Pose pose;
	protocols::farna::make_extended_coarse_pose( pose, fasta_sequence.sequence() );

	// native?
	pose::PoseOP native_pose;
	if ( option[ in::file::native ].user() ) {
		std::string native_pdb_file  = option[ in::file::native ];
		native_pose = new pose::Pose;
		import_pose::pose_from_pdb( *native_pose, *rsd_set_coarse, native_pdb_file );
	}

		//Constraints?
	if ( option[ cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
	}


	CoarseRNA_DeNovoProtocol de_novo_protocol( option[ out::nstruct ](), option[ cycles ](), option[ out::file::silent]() );
	de_novo_protocol.set_native_pose( native_pose );
	de_novo_protocol.set_temperature( option[ temperature]() );
	if( option[ score::weights ].user() ) de_novo_protocol.set_lores_scorefxn( option[ score::weights  ]() );
	if ( option[ params_file ].user() ) de_novo_protocol.set_rna_params_file( option[ params_file ]() );
	if ( option[data_file].user() )	de_novo_protocol.set_rna_data_file( option[ data_file ]() );
	de_novo_protocol.set_sim_anneal( !option[ no_sim_anneal ]() );
	de_novo_protocol.set_staged_constraints( option[ staged_constraints ]() );
	de_novo_protocol.set_close_loops( option[ close_loops ]() );
	de_novo_protocol.set_choose_best_solution( option[ choose_best_solution ]() );
	de_novo_protocol.set_force_ideal_chainbreak( option[ force_ideal_chainbreak ]() );
	de_novo_protocol.set_check_pairing_dists( option[ check_pairing_dists ]() );
	de_novo_protocol.set_add_base_pair_constraints( !option[skip_base_pair_constraints ]() );
	de_novo_protocol.set_dump_pdb( option[ dump ]() );
	de_novo_protocol.set_input_res( option[ input_res ]() );
	de_novo_protocol.set_freeze_domains( option[ freeze_domains ]() );

	if ( option[ in::file::silent ].user()	) de_novo_protocol.set_chunk_silent_files( option[ in::file::silent ]() );

	// visualize it.
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	de_novo_protocol.apply( pose );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	if ( option[ frag_test ]() ) {
		coarse_frag_test();
	}	else if ( option[ icoor ]() ) {
		icoor_test();
	} else if ( option[ pdbstats ]() ) {
		pdbstats_test();
	} else if ( option[ convert ]() ) {
		convert_to_coarse_test();
	} else if ( option[ output_coords ]() ) {
		output_minipose_coords_test();
	} else if ( option[ jump_database ]() ) {
		create_bp_jump_database_test();
	} else if ( option[ close_loop_test ]() ) {
		coarse_close_loop_test();
	} else if ( option[ rb_test ]() ) {
		coarse_rb_test();
	} else if ( option[ enumerate ]() ) {
		enumerate_map_test();
	} else if ( option[ little_motif ]() ) {
		pdb_little_motif_test();
	} else if ( option[ tar_motif ]() ) {
		tar_motif_test();
	} else if ( option[ mismatch ]() ) {
		mismatch_test();
	} else if ( option[ coarse_to_full ]() ) {
	  coarse_to_full_test();
	} else {
		coarse_rna_denovo_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}



///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;

	utility::vector1< Size > blank_size_vector;

	//Uh, options?
	NEW_OPT( frag_test, "blah",false);
	NEW_OPT( icoor, "blah",false);
	NEW_OPT( pdbstats, "blah",false);
	NEW_OPT( convert, "blah",false);
	NEW_OPT( output_coords, "blah",false);
	NEW_OPT( jump_database, "blah",false);
	NEW_OPT( no_sim_anneal, "blah",false);
	NEW_OPT( staged_constraints, "blah",false);
	NEW_OPT( close_loop_test, "blah",false);
	NEW_OPT( close_loops, "blah",false);
	NEW_OPT( choose_best_solution, "blah",false);
	NEW_OPT( force_ideal_chainbreak, "blah",false);
	NEW_OPT( check_pairing_dists, "blah",false);
	NEW_OPT( skip_base_pair_constraints, "blah",false);
	NEW_OPT( dump, "blah",false);
	NEW_OPT( rb_test, "blah",false);
	NEW_OPT( enumerate, "blah",false);
	NEW_OPT( sample_angles, "blah",false);
	NEW_OPT( output_vdw_pose, "blah",false);
	NEW_OPT( little_motif, "blah",false);
	NEW_OPT( tar_motif, "blah",false);
	NEW_OPT( mismatch, "blah",false);
	NEW_OPT( freeze_domains, "blah",false);
	NEW_OPT( coarse_to_full, "blah",false);
	NEW_OPT( cycles, "Default number of Monte Carlo cycles", 10000 );
	NEW_OPT( nbulge, "No. bulge residues", 3 );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	NEW_OPT( data_file, "Input file for RNA exposure data", "" );
	NEW_OPT( temperature, "Monte Carlo temperature", 5.0 );
	NEW_OPT( bin_width, "Monte Carlo temperature", 2.0 );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
	NEW_OPT( input_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector ); //I am here.

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
