¯// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/farna/RNA_StructureParameters.fwd.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_ChunkLibrary.hh>
#include <protocols/farna/RNA_ChunkLibrary.fwd.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/modeler/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/farna/RNA_LoopCloser.hh>
#include <protocols/farna/RNA_Minimizer.hh>

//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/util/basic.hh>
#include <core/io/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap

#include <core/scoring/Energies.hh> //for EnergyMap


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, build_helix_test )
OPT_KEY( Boolean, build_helix_precompute )
OPT_KEY( Boolean, all_combinations )
OPT_KEY( Boolean, minimize )
OPT_KEY( String,  params_file )
OPT_KEY( String,  seq )
OPT_KEY( String,  cst_file )
OPT_KEY( Integer, job_number )
OPT_KEY( Integer, total_jobs )

using core::pose::ResMap;


/////////////////////////////////////////////////////////////////////////////////
void
rna_assemble_test()
{
 	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace protocols::farna;
	using namespace core::kinematics;
	using namespace core::scoring::constraints;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

	core::pose::Pose pose;

	/////////////////////////////////////////////////////////////////////////////
	// Initialize an RNA with extended sequence + base pairs.
	//  Perhaps use RNA_StructureParameters class.
	//Prepare starting structure from scratch --> read from fasta.

	std::string seq_in;
	if ( option[ in::file::fasta ].user() ) {
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		seq_in = fasta_sequence->sequence();
	} else if ( option[ seq ].user() ){
		seq_in = option[ seq ]();
	} else {
		utility_exit_with_message( "Need to specify -fasta <fasta file> or -seq <sequence>." );
	}

	make_pose_from_sequence( pose,	seq_in,	*rsd_set );

	RNA_StructureParametersOP rna_structure_parameters( new RNA_StructureParameters );
	std::string jump_library_file( io::database::full_name("chemical/modeler/1jj2_RNA_jump_library.dat" ) );
	std::string rna_params_file( 	option[ params_file ] );

	rna_structure_parameters->initialize( pose, rna_params_file, jump_library_file, false );

	// Fill in helical segments with A-form jumps and torsion angles.
	// This will be new!!! Can do it last.

	//	rna_structure_parameters->setup_jumps( pose, true /*initialize_jumps*/ );
	rna_structure_parameters->setup_fold_tree_and_jumps_and_variants( pose );

	pose.dump_pdb( "init.pdb");

	if ( option[ cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	utility::vector1< std::string > const silent_files(  option[ in::file::silent ]() );
	RNA_ChunkLibrary rna_chunk_library( silent_files, pose.sequence(), rna_structure_parameters->connections() );

	for ( Size n = 1; n <= rna_chunk_library.num_chunk_sets(); n++ ) {
		//		std::cout << "Inserting chunk from library " << n << " into pose. " << std::endl;
		rna_chunk_library.insert_chunk_into_pose( pose, n, 1 );
		pose.dump_pdb( "blah"+string_of(n)+".pdb" );
	}
	// output.

	pose.dump_pdb( "chimera.pdb" );


	protocols::farna::RNA_LoopCloser rna_loop_closer;

	//	ScoreFunctionOP const lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	rna_loop_closer.close_loops_carefully( pose, rna_structure_parameters->connections() );

	pose.dump_pdb( "closed.pdb" );

	/////////////////////////////////////////////////
	if ( option[ out::file::silent ].user() ) {
		std::string tag( "blah" );
		std::string const silent_file = option[ out::file::silent  ]();
		SilentFileData silent_file_data;
		BinarySilentStruct s( pose, tag );
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
	}


}


/////////////////////////////////////////////////
std::string
get_tag( utility::vector1< Size > const & chunk_numbers ) {
	std::string tag( "S" );

	for ( Size i = 1; i <= chunk_numbers.size(); i++ ) {
		tag += "_";
		tag += string_of( chunk_numbers[i] );
	}

	//tag += ".pdb";
	return tag;
}

/////////////////////////////////////////////////
void
score_and_minimize( pose::Pose & pose, pose::Pose const & native_pose,
										std::string const & tag,
										std::string const & silent_file,
										protocols::farna::RNA_ChunkLibrary const & //rna_chunk_library
										)
{

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	static ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	lores_scorefxn->set_weight( atom_pair_constraint, 1.0 );

	(*lores_scorefxn)( pose );
	Real const atom_pair_constraint_score = (pose.energies()).total_energies()[ atom_pair_constraint ];
	//	std::cout << tag << " --> " << atom_pair_constraint_score << std::endl;

	static SilentFileData silent_file_data;
	static protocols::farna::RNA_Minimizer rna_minimizer;
	rna_minimizer.skip_o2prime_trials( true );
	//	rna_minimizer.set_allow_insert( rna_chunk_library.allow_insert() );

	if ( option[ minimize ]() ) {
		if ( atom_pair_constraint_score < 0.0 ) { /*arbitrary cutoff!!!*/

			pose::Pose pose_save( pose );

			rna_minimizer.apply( pose );

			BinarySilentStruct s( pose, tag );
			Real const rmsd = all_atom_rmsd( native_pose, pose );
			s.add_energy( "rms", rmsd );
			silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

			pose = pose_save;
		}
	} else {
		BinarySilentStruct s( pose, tag );
		Real const rmsd = all_atom_rmsd( native_pose, pose );
		s.add_energy( "rms", rmsd );
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}


}

/////////////////////////////////////////////////
void
insert_chunk( pose::Pose & pose, protocols::farna::RNA_ChunkLibrary const & rna_chunk_library,
							Size const chunk_set_index,
							utility::vector1 < Size > const & chunk_numbers_in,
							std::string const & silent_file,
							pose::Pose const & native_pose )
{

	using namespace core::options;
	using namespace core::options::OptionKeys;


	static Size const total_jobs_user = option[ total_jobs ];
	static Size const job_number_user = option[ job_number ];

	for ( Size n = 1; n <= rna_chunk_library.num_chunks( chunk_set_index ); n++ ) {

		if (chunk_set_index == 1 && ( n % total_jobs_user != job_number_user ) ) continue;

		//		std::cout << "Inserting chunk " << n << " from library " << chunk_set_index << " into pose. " << std::endl;
		rna_chunk_library.insert_chunk_into_pose( pose, chunk_set_index, n );

		utility::vector1 < Size > chunk_numbers( chunk_numbers_in );
		chunk_numbers.push_back( n );

		if ( chunk_set_index < rna_chunk_library.num_chunk_sets() ) {
			insert_chunk( pose, rna_chunk_library, chunk_set_index + 1, chunk_numbers, silent_file, native_pose );
		} else {

			std::string const tag = get_tag( chunk_numbers );
			//			std::cout << "Dumping: " << tag << std::endl;
			//			pose.dump_pdb( tag+".pdb" );

			score_and_minimize( pose, native_pose, tag, silent_file, rna_chunk_library );

		}
	}

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
void
rna_assemble_all_combinations_test()
{
 	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::io::silent;
	using namespace protocols::farna;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

	core::pose::Pose pose;

	/////////////////////////////////////////////////////////////////////////////
	// Initialize an RNA with extended sequence + base pairs.
	//  Perhaps use RNA_StructureParameters class.
	//Prepare starting structure from scratch --> read from fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	make_pose_from_sequence( pose,	fasta_sequence->sequence(),	*rsd_set );

	RNA_StructureParametersOP rna_structure_parameters( new RNA_StructureParameters );
	std::string jump_library_file( io::database::full_name("chemical/modeler/1jj2_RNA_jump_library.dat" ) );
	std::string rna_params_file( 	option[ params_file ] );

	rna_structure_parameters->initialize( pose, rna_params_file, jump_library_file, false );

	// Fill in helical segments with A-form jumps and torsion angles.
	// This will be new!!! Can do it last.
	//	rna_structure_parameters->setup_jumps( pose, true /*initialize_jumps*/ );
	rna_structure_parameters->setup_fold_tree_and_jumps_and_variants( pose );

	pose.dump_pdb( "init.pdb");

	utility::vector1< std::string > const silent_files(  option[ in::file::silent ]() );
	RNA_ChunkLibrary rna_chunk_library( silent_files, pose.sequence(), rna_structure_parameters->connections() );

	pose::Pose native_pose;
	std::string native_pdb_file  = option[ in::file::native ];
	io::pdb::pose_from_pdb( native_pose, *rsd_set, native_pdb_file );

	using namespace core::scoring::constraints;
	if ( option[ cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
	}

	std::string const silent_file = option[ out::file::silent  ]();

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	// RECURSIVELY TRY ALL COMBINATIONS!
	Size chunk_set_index( 1 );
	utility::vector1< Size > chunk_numbers;
	insert_chunk( pose, rna_chunk_library, chunk_set_index, chunk_numbers, silent_file, native_pose );

	pose.dump_pdb( "chimera1.pdb" );

}

/////////////////////////////////////////////////
void
rna_build_helix_test_OLD(){

 	using namespace core::scoring;
 	using namespace core::chemical;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace protocols::farna;

	std::string full_sequence;
	if ( option[ in::file::fasta ].user() ) {
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		full_sequence = fasta_sequence->sequence();
	} else if ( option[ seq ].user() ){
		full_sequence = option[ seq ]();
	} else {
		utility_exit_with_message( "Need to specify -fasta <fasta file> or -seq <sequence>." );
	}

	std::string silent_file = "";
	bool output_silent( false );
	if ( option[ out::file::silent ].user() ) {
		silent_file = option[ out::file::silent  ]();
		output_silent = true;
	}
	SilentFileData silent_file_data;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	make_pose_from_sequence( pose, full_sequence,	*rsd_set );

	RNA_HelixAssembler rna_helix_assembler;
	//rna_helix_assembler.random_perturbation( true );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Size const nstruct = option[ out::nstruct ];

	for (Size n = 1; n <= nstruct; n++ ) {
		rna_helix_assembler.apply( pose, full_sequence );

		std::string const tag( "S_"+lead_zero_string_of(n, 3) );
		if ( output_silent ) {
			BinarySilentStruct s( pose, tag );
			silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
		}
		pose.dump_pdb( full_sequence+".pdb" );
	}

}



/////////////////////////////////////////////////
void
rna_build_helix_test_precompute(){

 	using namespace core::scoring;
 	using namespace core::chemical;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::conformation;
	using namespace core::io::silent;
	using namespace core::io::pose_stream;
	using namespace protocols::farna;

	std::string full_sequence;
	if ( option[ in::file::fasta ].user() ) {
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		full_sequence = fasta_sequence->sequence();
	} else if ( option[ seq ].user() ){
		full_sequence = option[ seq ]();
	} else {
		utility_exit_with_message( "Need to specify -fasta <fasta file> or -seq <sequence>." );
	}

	std::string silent_file = "";
	bool output_silent( false );
	if ( option[ out::file::silent ].user() ) {
		silent_file = option[ out::file::silent  ]();
		output_silent = true;
	}
	SilentFileData silent_file_data;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	//	make_pose_from_sequence( pose, full_sequence,	*rsd_set );


	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// Move following to its own class when ready.
	///////////////////////////////////////////////////////////

	// six kinds of base pairs, so 36 kinds of 'base_pair_to_base_pair' files.
	// read them in from database.
	std::map< std::string, PoseOP > base_doublet_pose, base_pair_to_base_pair_pose;
	utility::vector1< std::string > base_doublets;
	base_doublets.push_back( "a_u" );
	base_doublets.push_back( "c_g" );
	base_doublets.push_back( "g_c" );
	base_doublets.push_back( "g_u" );
	base_doublets.push_back( "u_a" );
	base_doublets.push_back( "u_g" );

	for ( Size i = 1; i <= base_doublets.size(); i++ ){
		PoseOP pose = new Pose;
		std::string const outfile  = base_doublets[ i ] + "_anti_anti_north_north_antiparallel_00001.out";
		SilentFilePoseInputStream stream( outfile );
		stream.fill_pose( *pose, *rsd_set );
		base_doublet_pose[ pose->sequence() ] = pose;
	}

	for ( Size i = 1; i <= base_doublets.size(); i++ ){
		for ( Size j = 1; j <= base_doublets.size(); j++ ){
			PoseOP pose = new Pose;
			std::string const outfile  = base_doublets[ i ] + "_anti_anti_north_north_antiparallel_00001" + "_TO_" +
				base_doublets[ j ] +"_anti_anti_north_north_antiparallel_00001.out";
			SilentFilePoseInputStream stream( outfile );
			stream.fill_pose( *pose, *rsd_set );
			std::cout << "got it: " << pose->sequence() << std::endl;
			base_pair_to_base_pair_pose[ pose->sequence() ] = pose;
		}
	}

	Size nres = full_sequence.size();

	// starting doublet
	std::string sequence_start = full_sequence.substr(0,1) + full_sequence.substr(nres-1,1);
	make_pose_from_sequence( pose, sequence_start,	*rsd_set );
	FoldTree f( 2 );
	f.new_jump(1,2,1);
	f.set_jump_atoms( 1,
										core::chemical::rna::chi1_torsion_atom( pose.residue(1) ),
										core::chemical::rna::chi1_torsion_atom( pose.residue(2) ) );
	pose.fold_tree( f );
	PoseOP reference_pose_start = base_doublet_pose[ sequence_start ];

	ResMap res_map;
	res_map[1] = 1;
	res_map[2] = 2;
	std::cout << "SEQUENCE! " << pose.sequence() << std::endl;
	std::cout << "SEQ_REF!  " << reference_pose_start->sequence() << std::endl;
	std::cout << "BLAH!" << std::endl;
	copy_dofs_match_atom_names( pose, *reference_pose_start, res_map );
	std::cout << "PAST COPY_DOFS" << std::endl;

	for ( Size n = 2; n <= nres/2; n++ ){
		/////////////////////////////////////
		//this will become 'add base pair'
		/////////////////////////////////////
		char const seq1 = full_sequence[n-1];
		char const seq2 = full_sequence[nres-n];

		ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( seq1 ) )[1] ) ) );
		pose.append_polymer_residue_after_seqpos(   *rsd1, n - 1, true /*build_ideal_geometry*/ );

		ResidueOP rsd2( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( seq2 ) )[1] ) ) );
		pose.insert_residue_by_jump( *rsd2, n+1, n,
																 core::chemical::rna::chi1_torsion_atom( *rsd2 ),
																 core::chemical::rna::chi1_torsion_atom( *rsd1 ) );


		std::string base_pair_to_base_pair_seq = full_sequence.substr( n-2, 1)	+ full_sequence.substr(n-1,1)
			+ full_sequence.substr(nres-n,1) + full_sequence.substr(nres-n+1,1);
		PoseOP reference_pose = base_pair_to_base_pair_pose[ base_pair_to_base_pair_seq ];

		res_map.clear();
		res_map[ n-1 ] = 1;
		res_map[ n   ] = 2;
		res_map[ n+1 ] = 3;
		res_map[ n+2 ] = 4;

		std::cout << "SEQUENCE! " << std::endl;
		std::cout << "SEQ!    " << pose.sequence() << std::endl;
		std::cout << "BPSEQ!  " << base_pair_to_base_pair_seq << std::endl;
		std::cout << "REF_SEQ " << reference_pose->sequence() << std::endl;

		copy_dofs_match_atom_names( pose, *reference_pose, res_map );
	}

	pose.dump_pdb( "BLAH.pdb" );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ build_helix_test ] ){
		rna_build_helix_test_OLD();
	} else if ( option[ build_helix_precompute ] ){
		rna_build_helix_test_precompute();
	}	else if ( option[ all_combinations ] ){
			rna_assemble_all_combinations_test();
	} else {
		rna_assemble_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( build_helix_test, "build_helix_test", false );
	NEW_OPT( build_helix_precompute, "build_helix_test_precompute", false );
	NEW_OPT( minimize, "minimize", false );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	NEW_OPT( seq, "Input sequence", "" );
	NEW_OPT( all_combinations, "Systematically try all combinations", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
	NEW_OPT( job_number, "Job number", 0 );
	NEW_OPT( total_jobs, "Total jobs", 1 );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
