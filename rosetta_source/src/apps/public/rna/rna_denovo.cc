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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/RNA_DeNovoProtocol.hh>
#include <protocols/rna/RNA_DataReader.hh>
#include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>


// C++ headers
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/io/mpistream.hh>



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
OPT_KEY( Boolean, output_lores_silent_file )
OPT_KEY( Boolean, ignore_secstruct )
OPT_KEY( Boolean, filter_lores_base_pairs )
OPT_KEY( String, lores_scorefxn )
OPT_KEY( Boolean, vary_geometry )
OPT_KEY( Boolean, heat )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, staged_constraints )
OPT_KEY( Real, temperature )
OPT_KEY( Integer, cycles )
OPT_KEY( Real, jump_change_frequency )
OPT_KEY( String,  vall_torsions )
OPT_KEY( String,  jump_library_file )
OPT_KEY( String,  params_file )
OPT_KEY( String,  data_file )
OPT_KEY( String,  cst_file )


///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::rna;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Some initialization
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

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
		ensure_phosphate_nomenclature_matches_mini( native_pose );
		//	dump_pdb( native_pose, "native.pdb");
	}


	//Prepare starting structure from scratch --> read from fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( in_path + fasta_file )[1];
	core::pose::make_pose_from_sequence( extended_pose,	fasta_sequence->sequence(),	*rsd_set );

	//	dump_pdb( extended_pose, "extended.pdb");

	if (native_exists) std::cout << "Check it! NATIVE " << native_pose.sequence() << std::endl;
	std::cout << "Check it! EXTEND " << extended_pose.sequence() << std::endl;

	//Score these suckers.
	pose::Pose pose;
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// The good stuff
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	//Read in pose with ideal bond lengths and angles if it exists.
	pose = extended_pose;

	//	set_ideal_geometry( pose, extended_pose, rsd_set ); //by default, does nothing.
	bool input_pose( false );
	if ( option[in::file::s].user() ) {
		std::string const pdb_file = option[in::file::s][1];
		core::import_pose::pose_from_pdb( pose, *rsd_set, in_path + pdb_file );
		ensure_phosphate_nomenclature_matches_mini( pose );
		input_pose = true;
	}

	//	dump_pdb( pose, "start.pdb" );
	//	(*scorefxn)(pose);
	//	scorefxn->show( std::cout,pose );

	//Read in Torsion Library. Here go ahead and use my personal fragment class
	// because other instances of fragments in mini-rosetta (e.g., loop-modeling or ab initio)
	// are protein-specific!
	Size const nstruct = option[ out::nstruct ];
	Size const monte_carlo_cycles = option[ cycles ];

	std::string const silent_file = option[ out::file::silent  ]();
	bool const heat_structure( !input_pose || option[ heat ] );
	bool const minimize_structure = option[ minimize_rna ];
	bool const relax_structure = option[ relax_rna ];

	protocols::rna::RNA_DeNovoProtocol rna_de_novo_protocol( nstruct, monte_carlo_cycles,
																													 silent_file,
																													 heat_structure,
																													 minimize_structure,
																													 relax_structure );

	if (native_exists) rna_de_novo_protocol.set_native_pose( native_pose_OP );

	if ( option[ jump_library_file ].user() )	rna_de_novo_protocol.set_jump_library_file( in_path + option[ jump_library_file] );
 	if ( option[ vall_torsions ].user() )	rna_de_novo_protocol.set_vall_torsions_file( in_path + option[ vall_torsions ] );
	if ( option[params_file].user() )	rna_de_novo_protocol.set_rna_params_file( in_path + option[ params_file ] );
	if ( option[data_file].user() )	rna_de_novo_protocol.set_rna_data_file( in_path + option[ data_file ] );
	if ( option[lores_scorefxn].user() )	rna_de_novo_protocol.set_lores_scorefxn( option[ lores_scorefxn ] );

	rna_de_novo_protocol.ignore_secstruct( option[ ignore_secstruct ] );
	rna_de_novo_protocol.jump_change_frequency( option[ jump_change_frequency ] );
	rna_de_novo_protocol.close_loops( option[ close_loops] );
	rna_de_novo_protocol.output_lores_silent_file( option[ output_lores_silent_file ] );
	rna_de_novo_protocol.set_dump_pdb( option[ dump ] ) ;
	rna_de_novo_protocol.set_staged_constraints( option[ staged_constraints ] ) ;
	rna_de_novo_protocol.set_filter_lores_base_pairs(  option[ filter_lores_base_pairs] );
	rna_de_novo_protocol.set_vary_bond_geometry(  option[ vary_geometry ] );
	if ( option[ in::file::silent_struct_type ]() == "binary_rna" ) {
		rna_de_novo_protocol.set_binary_rna_output( true );
	}
	rna_de_novo_protocol.simple_rmsd_cutoff_relax( option[ simple_relax ] );
	if ( option[ in::file::silent ].user()	) {
		rna_de_novo_protocol.set_chunk_silent_files( option[ in::file::silent ]() );
	}

		//Constraints?
	if ( option[ cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[cst_file], new ConstraintSet, pose );
		pose.constraint_set( cst_set );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	rna_de_novo_protocol.apply( pose );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_denovo_test();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace basic::options;

	//Uh, options? MOVE THESE TO OPTIONS NAMESPACE INSIDE CORE/OPTIONS.
	NEW_OPT( minimize_rna, "Minimize RNA after fragment assembly",false );
	NEW_OPT( relax_rna, "Relax RNA after fragment assembly",false );
	NEW_OPT( simple_relax, "Relax by minimizing after any fragment insertion",false );
	NEW_OPT( ignore_secstruct, "Ignore sec struct in input file",false );
	NEW_OPT( lores_scorefxn, "Low resolution scorefunction weights file", "rna_lores.wts" );
	NEW_OPT( filter_lores_base_pairs, "Filter for models that satisfy structure parameters",false );
	NEW_OPT( vary_geometry, "Let bond lengths and angles vary from ideal",false );
	NEW_OPT( cycles, "Default number of Monte Carlo cycles", 10000 );
	NEW_OPT( vall_torsions, "Torsions file?", "1jj2.torsions" );
	NEW_OPT( temperature, "temperature", 0.3 );
	NEW_OPT( jump_change_frequency, "jump change frequency", 0.1 );
	NEW_OPT( close_loops, "close loops during frag insertion and jump mover", false );
	NEW_OPT( output_lores_silent_file, "output lores stuff", false );
	NEW_OPT( heat, "Heat (random frag insertions)", false );
	NEW_OPT( dump, "Dump pdb", false );
	NEW_OPT( staged_constraints, "Apply constraints in stages depending on sequence separation", false );
	NEW_OPT( jump_library_file, "Input file for jumps", "1jj2_RNA_jump_library.dat" );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	NEW_OPT( data_file, "Input file for RNA exposure data", "" );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

}
