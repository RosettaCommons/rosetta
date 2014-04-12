// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @authors Rhiju Das

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/util.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>
//#include <protocols/moves/PyMolMover.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, minimize_rna )
OPT_KEY( Boolean, relax_rna )
OPT_KEY( Boolean, simple_relax )
OPT_KEY( Boolean, close_loops )
OPT_KEY( Boolean, close_loops_after_each_move )
OPT_KEY( Boolean, output_lores_silent_file )
OPT_KEY( Boolean, binary_output )
OPT_KEY( Boolean, ignore_secstruct )
OPT_KEY( Boolean, filter_lores_base_pairs )
OPT_KEY( Boolean, filter_lores_base_pairs_early )
OPT_KEY( Boolean, filter_chain_closure )
OPT_KEY( Boolean, use_1jj2_torsions )
OPT_KEY( String, lores_scorefxn )
OPT_KEY( Boolean, heat )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, staged_constraints )
OPT_KEY( Real, temperature )
OPT_KEY( Real, rna_lores_linear_chainbreak_weight )
OPT_KEY( Real, rna_lores_chainbreak_weight )
OPT_KEY( Integer, cycles )
OPT_KEY( Real, jump_change_frequency )
OPT_KEY( Real, suppress_bp_constraint )
OPT_KEY( String,  refine_silent_file )
OPT_KEY( String,  jump_library_file )
OPT_KEY( String,  params_file )
OPT_KEY( String,  data_file )
OPT_KEY( String,  cst_file )
OPT_KEY( IntegerVector, chunk_res ) // deprecated
OPT_KEY( Boolean, allow_bulge )
OPT_KEY( Boolean, allow_consecutive_bulges )
OPT_KEY( IntegerVector, allowed_bulge_res )
OPT_KEY( IntegerVector, extra_minimize_res )
OPT_KEY( Boolean, root_at_first_rigid_body )
OPT_KEY( Boolean, move_first_rigid_body )
OPT_KEY( Boolean, output_filters )
OPT_KEY( Boolean, autofilter )
OPT_KEY( Real, filter_chain_closure_distance )
OPT_KEY( Boolean, filter_chain_closure_halfway )
OPT_KEY( IntegerVector, output_res_num )
OPT_KEY( Boolean, refine_native )
OPT_KEY( Boolean, bps_moves )
OPT_KEY( Boolean, minimizer_use_coordinate_constraints )


///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::farna;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Some initialization
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::PoseOP native_pose_OP = new pose::Pose;
	pose::Pose & native_pose = *native_pose_OP;

	pose::Pose extended_pose;

	std::string const in_path = option[ in::path::path ]()[1];

	bool native_exists = false;
	//Read in native if it exists.
	if ( option[ in::file::native ].user() ) {
		native_exists = true;
		//Read in native if it exists.
		std::string native_pdb_file  = option[ in::file::native ];
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, in_path + native_pdb_file );
		// ensure_phosphate_nomenclature_matches_mini( native_pose );
		//dump_pdb( native_pose, "native.pdb");
	} else {
		runtime_assert( !option[refine_native]() );
	}


	//Prepare starting structure from scratch --> read from fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( in_path + fasta_file )[1];
	if ( option[refine_native]() ) {
		extended_pose = native_pose;
	} else {
		core::pose::make_pose_from_sequence( extended_pose,	fasta_sequence->sequence(),	*rsd_set );
	}
	set_output_res_num( extended_pose, option[ output_res_num ]() );

	// Silent file input for fine refinement
	utility::vector1<pose::PoseOP> refine_pose_list;
	if ( option[refine_silent_file]() != "" ) {
		core::import_pose::pose_stream::SilentFilePoseInputStream input( option[refine_silent_file]() );
		input.set_order_by_energy( true );
		while ( input.has_another_pose() ) {
			pose::PoseOP new_pose = new pose::Pose;
			input.fill_pose( *new_pose, *rsd_set );
			set_output_res_num( *new_pose, option[ output_res_num ]() );
			refine_pose_list.push_back( new_pose );
		}
	}

	//	dump_pdb( extended_pose, "extended.pdb");

	if (native_exists) std::cout << "Check it! NATIVE " << native_pose.sequence() << std::endl;
	std::cout << "Check it! EXTEND " << extended_pose.sequence() << std::endl;

	//Score these suckers.
	pose::Pose pose( extended_pose );
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// The good stuff
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	//Read in pose with ideal bond lengths and angles if it exists.

	//	set_ideal_geometry( pose, extended_pose, rsd_set ); //by default, does nothing.

	//Read in Torsion Library. Here go ahead and use my personal fragment class
	// because other instances of fragments in mini-rosetta (e.g., loop-modeling or ab initio)
	// are protein-specific!
	Size const nstruct = option[ out::nstruct ];
	std::string const silent_file = option[ out::file::silent  ]();
	bool heat_structure( true );
	if ( option[refine_native]() || option[refine_silent_file]() != "" ) heat_structure = false;
	bool const minimize_structure = option[ minimize_rna ];
	bool const relax_structure = option[ relax_rna ];
	bool const is_allow_bulge = option[ allow_bulge ];

  protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol( nstruct,
                                                           silent_file,
                                                           heat_structure,
                                                           minimize_structure,
                                                           relax_structure,
                                                           is_allow_bulge );

	if (native_exists) rna_de_novo_protocol.set_native_pose( native_pose_OP );
	if ( option[ cycles ].user() )	rna_de_novo_protocol.set_monte_carlo_cycles( option[ cycles ]() );
	if ( option[ jump_library_file ].user() )	rna_de_novo_protocol.set_jump_library_file( in_path + option[ jump_library_file] );
 	if ( option[ basic::options::OptionKeys::rna::vall_torsions ].user() )	{
		// check in database first
		std::string vall_torsions_file = basic::database::full_name("/sampling/rna/" + option[ basic::options::OptionKeys::rna::vall_torsions ]() );
		if (!utility::file::file_exists( vall_torsions_file ) && !utility::file::file_exists( vall_torsions_file + ".gz" ) )  vall_torsions_file = in_path + option[ basic::options::OptionKeys::rna::vall_torsions ]();
		rna_de_novo_protocol.set_vall_torsions_file( vall_torsions_file );
	}
	if ( option[ use_1jj2_torsions ]() ) rna_de_novo_protocol.set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	if ( option[params_file].user() )	rna_de_novo_protocol.set_rna_params_file( in_path + option[ params_file ] );
	if ( option[data_file].user() )	rna_de_novo_protocol.set_rna_data_file( in_path + option[ data_file ] );
	if ( option[lores_scorefxn].user() )	rna_de_novo_protocol.set_lores_scorefxn( option[ lores_scorefxn ] );
	if ( option[ rna_lores_linear_chainbreak_weight ].user() ) rna_de_novo_protocol.set_linear_chainbreak_weight( option[ rna_lores_linear_chainbreak_weight ]() );
	if ( option[ rna_lores_chainbreak_weight ].user() ) rna_de_novo_protocol.set_chainbreak_weight( option[ rna_lores_chainbreak_weight ]() );

	rna_de_novo_protocol.set_temperature( option[ temperature ] );
	rna_de_novo_protocol.ignore_secstruct( option[ ignore_secstruct ] );
	rna_de_novo_protocol.jump_change_frequency( option[ jump_change_frequency ] );
	rna_de_novo_protocol.set_close_loops( option[ close_loops] );
	rna_de_novo_protocol.set_close_loops_after_each_move( option[ close_loops_after_each_move ] );
	rna_de_novo_protocol.output_lores_silent_file( option[ output_lores_silent_file ] );
	rna_de_novo_protocol.set_dump_pdb( option[ dump ] ) ;
	rna_de_novo_protocol.set_staged_constraints( option[ staged_constraints ] ) ;
	rna_de_novo_protocol.set_filter_lores_base_pairs(  option[ filter_lores_base_pairs] );
	rna_de_novo_protocol.set_filter_lores_base_pairs_early(  option[ filter_lores_base_pairs_early] );
	rna_de_novo_protocol.set_suppress_bp_constraint(  option[ suppress_bp_constraint]() );
	rna_de_novo_protocol.set_filter_chain_closure(  option[ filter_chain_closure ]() );
	rna_de_novo_protocol.set_filter_chain_closure_distance(  option[ filter_chain_closure_distance ]() );
	rna_de_novo_protocol.set_filter_chain_closure_halfway(  option[ filter_chain_closure_halfway ]() );
	rna_de_novo_protocol.set_vary_bond_geometry(  option[ basic::options::OptionKeys::rna::vary_geometry ] );
	rna_de_novo_protocol.set_move_first_rigid_body(  option[ move_first_rigid_body ] );
	rna_de_novo_protocol.set_root_at_first_rigid_body(  option[ root_at_first_rigid_body ] );
	rna_de_novo_protocol.set_output_filters(  option[ output_filters ] );
	rna_de_novo_protocol.set_autofilter(  option[ autofilter ] );
	rna_de_novo_protocol.set_bps_moves(  option[ bps_moves ] );
	rna_de_novo_protocol.set_minimizer_use_coordinate_constraints( option[ minimizer_use_coordinate_constraints ]() );
	if ( option[ in::file::silent_struct_type ]() == "binary_rna"  || option[ binary_output ]() )	rna_de_novo_protocol.set_binary_rna_output( true );

	rna_de_novo_protocol.simple_rmsd_cutoff_relax( option[ simple_relax ] );
	if ( option[ in::file::s ].user() ) rna_de_novo_protocol.set_chunk_pdb_files( option[ in::file::s ]() );
	if ( option[ in::file::silent ].user() ) 	rna_de_novo_protocol.set_chunk_silent_files( option[ in::file::silent ]() );

	runtime_assert( ! (option[ chunk_res].user() && option[ in::file::input_res ].user() ) ); // let's deprecate chunk_res soon
	if ( option[ chunk_res ].user() ) {
		std::cout << "WARNING! WARNING! WARNING! -chunk_res will be deprecated soon. Use -input_res instead!" << std::endl;
		std::cerr << "WARNING! WARNING! WARNING! -chunk_res will be deprecated soon. Use -input_res instead!" << std::endl;
		rna_de_novo_protocol.set_input_res( option[ chunk_res ]() ) ;
	}

	if ( option[ in::file::input_res ].user() )  rna_de_novo_protocol.set_input_res( option[ in::file::input_res ]() ) ;

	rna_de_novo_protocol.set_allow_consecutive_bulges( option[ allow_consecutive_bulges ]() ) ;
	rna_de_novo_protocol.set_allowed_bulge_res( option[ allowed_bulge_res ]() ) ;
	rna_de_novo_protocol.set_extra_minimize_res( option[ extra_minimize_res ]() ) ;

	rna_de_novo_protocol.set_refine_pose_list( refine_pose_list );
	if ( option[refine_native]() || option[refine_silent_file]() != "" ) rna_de_novo_protocol.set_rounds( 1 );
	if ( option[refine_native]() ) rna_de_novo_protocol.set_refine_pose( true );

	//Constraints?
	if ( option[ cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
		for ( Size i = 1; i <= refine_pose_list.size(); ++i) {
			// ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, *refine_pose_list[i] );
			refine_pose_list[i]->constraint_set( cst_set );
		}
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );

	//	protocols::moves::AddPyMolObserver( pose, false, 0.01);

	rna_de_novo_protocol.apply( pose );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_denovo_test();
  protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
try {
	using namespace basic::options;


	std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence>  [ -native <native pdb file> ] " << std::endl;
	std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

	utility::vector1< Size > blank_size_vector;

	option.add_relevant( in::file::input_res );

	NEW_OPT( minimize_rna, "Minimize RNA after fragment assembly",false );
	NEW_OPT( relax_rna, "Relax RNA after fragment assembly",false );
	NEW_OPT( simple_relax, "Relax by minimizing after any fragment insertion",false );
	NEW_OPT( ignore_secstruct, "Ignore sec struct in input file",false );
	NEW_OPT( lores_scorefxn, "Low resolution scorefunction weights file", "farna/rna_lores.wts" );
	NEW_OPT( filter_lores_base_pairs, "Filter for models that satisfy structure parameters", true );
	NEW_OPT( filter_lores_base_pairs_early, "Filter for models that satisfy structure parameters at round 2 of 10", true );
	NEW_OPT( filter_chain_closure, "Filter for models that have closed chains after lores before minimize",true );
	NEW_OPT( filter_chain_closure_halfway, "Filter for models that have closed chains after lores before minimize at round 5 of 10", true );
	NEW_OPT( filter_chain_closure_distance, "Mean distance across 3 chainbreak atoms to filter models that have closed chains after lores", 6.0 );
	NEW_OPT( cycles, "Default number of Monte Carlo cycles", 0 ); // now default is set based on the number of moving residues.
	NEW_OPT( temperature, "temperature", 2.0 );
	NEW_OPT( jump_change_frequency, "jump change frequency", 0.1 );
	NEW_OPT( close_loops, "close loops after de novo protocol and again after minimization", true ); /* this should always be on. */
	NEW_OPT( close_loops_after_each_move, "close loops during frag insertion and jump mover -- can be expensive", false );
	NEW_OPT( output_lores_silent_file, "output lores stuff", false );
	NEW_OPT( heat, "Heat (random frag insertions)", false );
	NEW_OPT( dump, "Dump pdb", false );
	NEW_OPT( staged_constraints, "Apply constraints in stages depending on sequence separation", false );
	NEW_OPT( jump_library_file, "Input file for jumps", "sampling/rna/1jj2_RNA_jump_library.dat" );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	NEW_OPT( data_file, "Input file for RNA exposure data", "" );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );
	NEW_OPT( chunk_res, "Input residues for chunk libraries (specified by -in:file:silent or -s) ... use -input_res instead!", blank_size_vector );
	NEW_OPT( use_1jj2_torsions, "Use original (ribosome) fragments, 1JJ2", false );
	NEW_OPT( rna_lores_chainbreak_weight, "chainbreak weight for lo res sampling", 0.0 );
	NEW_OPT( rna_lores_linear_chainbreak_weight, "linear chainbreak weight for lo res sampling", 0.0 );
  NEW_OPT( allow_bulge , "Automatically virtualize residues that are not energetically stable", false );
  NEW_OPT( allowed_bulge_res, "Use with allow_bulge, allowable pos for virtualization", blank_size_vector  );
  NEW_OPT( extra_minimize_res, "Extra residues during minimize step", blank_size_vector  );
  NEW_OPT( allow_consecutive_bulges, "allow_consecutive_bulges", false );
  NEW_OPT( binary_output, "force output to binary rna silentstruct", false );
  NEW_OPT( move_first_rigid_body, "first_rigid_body is usually kept frozen, but might be useful to sample it.", false );
  NEW_OPT( root_at_first_rigid_body, "places coordinate system away from the usual last virtual residue and puts it on the first rigid body. useful if this rigidbody needs to be fixed, but other bodies need to move as if this one is moving. Use with -move_first_rigid_body. ", false );
  NEW_OPT( suppress_bp_constraint, "Factor by which to lower base pair constraint weight. ", 1.0 );
  NEW_OPT( output_filters, "output lores scores at early stage (round  2 of 10) and at end -- could be useable for early termination of unpromising early starts", false );
  NEW_OPT( autofilter, "Automatically skip output/minimize if lores score is worse than 20th percentile, updated on the fly.", true );
  NEW_OPT( output_res_num, "Numbering of residues in output PDB or silent file", blank_size_vector  );
  NEW_OPT( refine_silent_file, "Name of the silent file to be refined.", "" );
  NEW_OPT( refine_native, "Refine starting from the native pose", false );
  NEW_OPT( bps_moves, "Base pair step moves", false );
  NEW_OPT( minimizer_use_coordinate_constraints, "Use coordinate constraints for first round of minimizer", true );

	option.add_relevant( basic::options::OptionKeys::rna::vary_geometry );
	option.add_relevant( basic::options::OptionKeys::rna::vall_torsions );
	//	option.add_relevant( basic::options::OptionKeys::rna::jump_database );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}

