// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   apps/pilot/ronj/find_position_matches_using_inverse_rotamers.cc
/// @brief  Program which scans for backbone positions in a list of scaffolds with low rmsd to 
///			specific positions in a reference structure by building inverse rotamers
/// @author Ron Jacak (ron.jacak@gmail.com)


#include <basic/options/util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

/// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

// Protocol headers

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <cstdio>
//#include <cstdlib>

static basic::Tracer TR("apps.pilot.ronj.find_position_matches_using_inverse_rotamers");

//#define FILE_DEBUG 1

// application specific options
namespace find_position_matches_using_inverse_rotamers {

	basic::options::FileOptionKey const reference_structure( "find_position_matches_using_inverse_rotamers::reference_structure" );

	basic::options::StringVectorOptionKey const scaffold_epitope_residue_range( "find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range" );
	basic::options::StringVectorOptionKey const reference_structure_epitope_residue_range( "find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range" );
	basic::options::StringVectorOptionKey const reference_structure_positions_to_match( "find_position_matches_using_inverse_rotamers::reference_structure_positions_to_match" );
	basic::options::StringVectorOptionKey const reference_structure_glycan_residues( "find_position_matches_using_inverse_rotamers::reference_structure_glycan_residues" );

	//basic::options::IntegerOptionKey const num_rotamers( "find_position_matches_using_inverse_rotamers::num_rotamers" );
	basic::options::RealOptionKey const ir_rmsd( "find_position_matches_using_inverse_rotamers::ir_rmsd" );
	basic::options::IntegerOptionKey const number_glycan_clashes_allowed( "find_position_matches_using_inverse_rotamers::number_glycan_clashes_allowed" );

	basic::options::BooleanOptionKey const output_inverse_rotamers( "find_position_matches_using_inverse_rotamers::output_inverse_rotamers" );
	basic::options::BooleanOptionKey const output_aligned_scaffold( "find_position_matches_using_inverse_rotamers::output_aligned_scaffold" );
	basic::options::BooleanOptionKey const check_for_valid_chi3_chi4_dihedrals( "find_position_matches_using_inverse_rotamers::check_for_valid_chi3_chi4_dihedrals" );
	basic::options::BooleanOptionKey const assume_command_line_has_pose_numbering( "find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering" );


} // end namespace find_position_matches_using_inverse_rotamers


std::string usage_string;

///
/// @begin init_usage_prompt
///
/// @brief
/// the usage prompt that gets printed when the user doesn't enter all the required command line arguments
///
void init_usage_prompt( std::string exe ) {

	// place the prompt up here so that it gets updated easily; global this way, but that's ok
	std::stringstream usage_stream;
	usage_stream
			<< "Usage: " << exe
			<< "\n\t-database path/to/rosetta/database"
			<< "\n\t-s pdb                                                 read in input PDB file or silent file"
			<< "\n\t"
			<< "\n\treference_structure                                    PDB file which contains antibody and glycosylated epitope"
			<< "\n\tscaffold_epitope_residue_range                         epitope chain and residue start and stop ids in scaffold, PDB numbering (e.g. A/47 A/51)"
			<< "\n\treference_structure_epitope_residue_range              epitope chain residue start and stop ids in reference structure, PDB numbering"
			<< "\n\treference_structure_positions_to_match                 chain and residue ids of positions to build inverse rotamers for, PDB numbering (e.g. D/501)"
			<< "\n\treference_structure_glycan_residues                    chain and residue ids of glycan molecules in reference structure, PDB numbering"

			<< "\n\tir_rmsd                                                   rmsd limit between scaffold position and inverse rotamer (default: 0.5)"
			<< "\n\tnumber_glycan_clashes_allowed                          maximum amount of clash between glycans and rest of scaffold (default: 100 )"

			<< "\n\toutput_inverse_rotamers                                output PDB files for the inverse rotamers that pass all filters (default: false)"
			<< "\n\toutput_aligned_scaffold                                output PDB file for scaffold aligned to reference structure (default: false)"
			<< "\n\tcheck_for_valid_chi3_chi4_dihedrals                    check that the chi3 and chi4 angles of matching rotamers are within observed glycan distribution values (default: false)"
			<< "\n\tassume_command_line_has_pose_numbering                 instead of translating the input 'A/47' to pose numbering, assume command line arguments are already pose numbering (default: false)"

			<< "\n\n";
	usage_string = usage_stream.str();

}


///
/// @begin build_bb_independent_rotamers
///
/// @brief
/// function which returns a set of backbone-independent rotamers for the given ResidueType.
/// basically a copy of the function in core::pack::rotamer_set::bb_independent_rotamers.hh, but it samples more chi angles
/// thereby building more rotamers. 
///
void build_bb_independent_rotamers( core::chemical::ResidueTypeCOP target_res_restype, utility::vector1< core::conformation::ResidueOP > & rotamers ) {

	using namespace core;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::pack::rotamer_set;
	using namespace core::graph;
	
	conformation::Residue target_res( *target_res_restype, true ); // boolean is a required dummy argument for Residue constructor
	pose::Pose dummy_pose;
	dummy_pose.append_residue_by_jump( target_res, (core::Size) 0 );
	pose::add_lower_terminus_type_to_pose_residue( dummy_pose, 1 /* seqpos */); // probably critical so that the dunbrack library uses neutral phi
	pose::add_upper_terminus_type_to_pose_residue( dummy_pose, 1 /* seqpos */); // probably critical so that the dunbrack library uses neutral psi

	scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn( dummy_pose );

	pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( dummy_pose );
	task->initialize_from_command_line(); // initialize rotamer options from the command line. this allows the user to toggle how many rotamers are built without the need for recompiling.
	task->nonconst_residue_task( 1 ).restrict_to_repacking(); // restrict to repacking call ensures that we only build rotamers of the native amino acid type
	task->nonconst_residue_task( 1 ).or_include_current( false ); // can't "include_current" because there is no current or native residue

	// do extra sampling around chi angles. these can be activating via the command line so leave them commented out here.
	//task->nonconst_residue_task(1).or_ex1( true );
	//task->nonconst_residue_task(1).or_ex2( true );
	//task->nonconst_residue_task(1).or_ex3( true );
	//task->nonconst_residue_task(1).or_ex4( true );
	//task->nonconst_residue_task(1).or_ex1_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	//task->nonconst_residue_task(1).or_ex2_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	//task->nonconst_residue_task(1).or_ex3_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
	//task->nonconst_residue_task(1).or_ex4_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);

	graph::GraphOP dummy_png = core::pack::create_packer_graph( dummy_pose, dummy_sfxn, task );

	pack::rotamer_set::RotamerSetFactory rsf;
	pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( dummy_pose.residue( 1 ) ) );
	rotset->set_resid( 1 );
	bool use_neighbor_context = false;
	rotset->build_rotamers( dummy_pose, dummy_sfxn, *task, dummy_png, use_neighbor_context );

	// now when creating the rotamers, we have to make sure we don't sneak in the additional variant types
	for( core::Size i = 1; i <= rotset->num_rotamers(); ++i ){
		core::conformation::ResidueOP rot( target_res.clone() );
		for( core::Size j = 1; j <= target_res.nchi(); ++j ) {
			rot->set_chi( j, rotset->rotamer( i )->chi( j ) );
		}
		rotamers.push_back( rot );
	}

	return;
}


void tokenize_string( const std::string & str, utility::vector1< std::string > & tokens, const std::string & delimiters = " " ) {

	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

	// Find first "non-delimiter".
	std::string::size_type pos = str.find_first_of(delimiters, lastPos);

	while ( std::string::npos != pos || std::string::npos != lastPos ) {
		// Found a token, add it to the vector.
		tokens.push_back( str.substr( lastPos, pos - lastPos ) );

		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of( delimiters, pos );

		// Find next "non-delimiter"
		pos = str.find_first_of( delimiters, lastPos );
	}
}

void parse_position( std::string position, char & chain, core::Size & pdb_resnum, char & icode ) {

	// tokenize the strings into resnum and chain (and icode, if applicable)
	utility::vector1< std::string > parsed_tokens;

	// tokenize the string using the delimeter slash /
	tokenize_string( position, parsed_tokens, "/" );

	// use an istringstream to convert the string chain id to char chain id
	std::istringstream iss( parsed_tokens[ 1 ] );
	iss >> chain;
	
	std::string position_code = parsed_tokens[ 2 ];
	
	// check to see if an insertion code is present in the position_code string
	// if the string is made of all digits, no icode is present
	std::stringstream ss;

	if ( position_code.find_first_not_of("0123456789") == std::string::npos ) {
		icode = ' ';
		ss << position_code;
		ss >> pdb_resnum;

	} else {
		for ( std::string::iterator it = position_code.begin(); it < position_code.end(); ++it ) {
			if ( isdigit(*it) ) {
				ss << (*it);
			} else {
				icode = *it; // assumes that insertion code is only 1-letter, and that no digits come after it!
			}
		}
		ss >> pdb_resnum; // converts the ss buffer contents to a Size type
	}

	return;
}


///
/// @begin dump_rotamerset_pdb
///
/// @brief
/// function which writes out a PDB file of the rotamers in the passed in vector
/// this should really be a function of the RotamerSet class, but it's not. there is a similar function in the class RotamerSets, 
/// but I don't want the multiple MODELS in there.
/// only using this function to check that the inverse rotamers are aligned properly to the reference residue(s).
///
void dump_rotamerset_pdb( utility::vector1< core::conformation::ResidueCOP > const & rotamers, std::string const & filename ) {

	using core::Size;

	// open file
	std::ofstream out( filename.c_str() );

	Size atom_counter = 0;
	for ( Size ii = 1; ii <= rotamers.size(); ++ii ) {
		core::io::pdb::dump_pdb_residue( *(rotamers[ ii ]), atom_counter, out );
	}

	out.close();
}


int
main( int argc, char * argv [] ) {

	try {

	using namespace basic::options;
	using namespace core;

	// can't use the option remodel::help because options don't get init'd until after devel::init. use argc/argv instead.
	if ( argc == 1 || ( argc > 1 && strcmp(argv[ 1 ], "-h") == 0 ) ) {
		init_usage_prompt( argv[0] );
		std::cout << usage_string;
		exit(0);
	}

	// add application specific options to options system
	option.add( find_position_matches_using_inverse_rotamers::reference_structure, "reference_structure" );

	option.add( find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range, "scaffold_epitope_residue_range" );
	option.add( find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range, "reference_structure_epitope_residue_range" );
	option.add( find_position_matches_using_inverse_rotamers::reference_structure_positions_to_match, "reference_structure_positions_to_match" );
	option.add( find_position_matches_using_inverse_rotamers::reference_structure_glycan_residues, "reference_structure_glycan_residues" );

	//option.add( find_position_matches_using_inverse_rotamers::num_rotamers, "num_rotamers" ).def( 10 );
	option.add( find_position_matches_using_inverse_rotamers::ir_rmsd, "ir_rmsd" ).def( 0.5 );
	option.add( find_position_matches_using_inverse_rotamers::number_glycan_clashes_allowed, "number_glycan_clashes_allowed" ).def( 100 );

	option.add( find_position_matches_using_inverse_rotamers::output_inverse_rotamers, "output_inverse_rotamers" ).def( false );
	option.add( find_position_matches_using_inverse_rotamers::output_aligned_scaffold, "output_aligned_scaffold" ).def( false );
	option.add( find_position_matches_using_inverse_rotamers::check_for_valid_chi3_chi4_dihedrals, "check_for_valid_chi3_chi4_dihedrals" ).def( false );

	option.add( find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering, "assume_command_line_has_pose_numbering" ).def( false );

	
	// initialize rosetta
	devel::init( argc, argv );

	// open an output file
	FILE* fout;
	fout = fopen( "single_matches.txt", "w" );

	// add a header line to the output file to let the user know what the fields are
	fprintf( fout, "chain-resi-icode\tconformation\tfilename\taligned_scaffold\tno.residues\tchain-resi-icode\tRMSD\tno.clashes\tchi1\tchi2\tchi3\tchi4\tdotp_angle\tdihedral\n" );
	fflush( fout );

	// the reference structure contains the antibody and the glycosylated epitope. 
	// read it in up here, outside of the loop of all scaffolds to match
	pose::Pose reference_structure;
	import_pose::pose_from_pdb( reference_structure, option[ find_position_matches_using_inverse_rotamers::reference_structure ]() );

	// need to convert the StringVectorOptionKey option to start and stop residues, because that's what we'll align on later.
	// we'll assume that the user specified the range as start to stop.
	if ( !option[ find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range ].user() ) {
		std::cerr << "ERROR: Epitope residue range on reference structure not specified." << std::endl;
		utility_exit();
	} else if ( option[ find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range ]().size() > 2 ) {
		std::cerr << "ERROR: Invalid epitope residue range on reference structure specified." << std::endl;
		utility_exit();
	}

	std::string reference_structure_epitope_residue_start_pdb = option[ find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range ]()[ 1 ];
	std::string reference_structure_epitope_residue_stop_pdb = option[ find_position_matches_using_inverse_rotamers::reference_structure_epitope_residue_range ]()[ 2 ];

	TR << "reference structure epitope residue range (PDB numbering): " << reference_structure_epitope_residue_start_pdb << "-" << reference_structure_epitope_residue_stop_pdb << std::endl;

	Size reference_structure_epitope_residue_start;
	Size reference_structure_epitope_residue_stop;

	if ( !option[ find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering ] ) {
		// convert PDB numbering to pose numbering, because we use the pose numbering later
		Size pdb_resnum; 
		char chain, icode;
		parse_position( reference_structure_epitope_residue_start_pdb, chain, pdb_resnum, icode );
		reference_structure_epitope_residue_start = (reference_structure.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

		parse_position( reference_structure_epitope_residue_stop_pdb, chain, pdb_resnum, icode );
		reference_structure_epitope_residue_stop = (reference_structure.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

	} else {
		std::istringstream iss( reference_structure_epitope_residue_start_pdb );
		iss >> reference_structure_epitope_residue_start;
		
		iss.clear();
		iss.str( reference_structure_epitope_residue_stop_pdb );
		iss >> reference_structure_epitope_residue_stop;
	}

	TR << "reference structure epitope residue range (pose numbering): " << reference_structure_epitope_residue_start << "-" << reference_structure_epitope_residue_stop << std::endl;
		
	// make sure the user specified asparagine positions; these are the positions we are trying to match for in the scaffolds
	if ( !option[ find_position_matches_using_inverse_rotamers::reference_structure_positions_to_match ].user() ) {
		std::cerr << "ERROR: No positions on reference structure to match specified." << std::endl;
		utility_exit();
	}
	utility::vector1< std::string > positions_to_match_pdb = option[ find_position_matches_using_inverse_rotamers::reference_structure_positions_to_match ];
	utility::vector1< Size > positions_to_match;
	for ( Size ii = 1; ii <= positions_to_match_pdb.size(); ++ii ) {

		if ( !option[ find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering ] ) {
			Size pdb_resnum;
			char chain, icode;
			parse_position( positions_to_match_pdb[ ii ], chain, pdb_resnum, icode );
			positions_to_match.push_back( (reference_structure.pdb_info())->pdb2pose( chain, pdb_resnum, icode ) );
		} else {
			std::istringstream iss( positions_to_match_pdb[ ii ] );
			Size position;
			iss >> position;
			positions_to_match.push_back( position );
			iss.clear();
		}
	}

	TR << "reference structure positions to match: ";
	for ( Size ii = 1; ii <= positions_to_match.size(); ++ii ) {
		TR << positions_to_match[ ii ];
		if ( ii != positions_to_match.size() ) TR << ", ";
	}
	TR << std::endl;

	// make sure the user specified glycan positions, too
	if ( !option[ find_position_matches_using_inverse_rotamers::reference_structure_glycan_residues ].user() ) {
		std::cerr << "ERROR: No glycan positions on reference structure specified." << std::endl;
		utility_exit();
	}
	utility::vector1< std::string > reference_structure_glycan_residues_pdb = option[ find_position_matches_using_inverse_rotamers::reference_structure_glycan_residues ];
	utility::vector1< Size > reference_structure_glycan_residues;
	for ( Size ii = 1; ii <= reference_structure_glycan_residues_pdb.size(); ++ii ) {
		if ( !option[ find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering ] ) {
			Size pdb_resnum;
			char chain, icode;
			parse_position( reference_structure_glycan_residues_pdb[ ii ], chain, pdb_resnum, icode );
			reference_structure_glycan_residues.push_back( (reference_structure.pdb_info())->pdb2pose( chain, pdb_resnum, icode ) );
		} else {
			std::istringstream iss( reference_structure_glycan_residues_pdb[ ii ] );
			Size position;
			iss >> position;
			reference_structure_glycan_residues.push_back( position );
			iss.clear();
		}
	}

	TR << "reference structure glycan residues: ";
	for ( Size ii = 1; ii <= reference_structure_glycan_residues.size(); ++ii ) {
		TR << reference_structure_glycan_residues[ ii ];
		if ( ii != reference_structure_glycan_residues.size() ) TR << ", ";
	}
	TR << std::endl;

	// save it to a variable to keep from having to check the option over and over
	Size number_glycan_clashes_allowed = option[ find_position_matches_using_inverse_rotamers::number_glycan_clashes_allowed ];

	//
	// concatenate -s and -l flags together to get total list of PDB files
	// The called function will die with a useful error message if neither -s or -l is specified.
	// Check to make sure all of the files exist here, too.
	//
	utility::vector1< std::string > pdb_file_names = basic::options::start_files();
	utility::vector1< std::string >::iterator iter_input_pdb_filename, iter_last_pdb;

	// now, for each "hybrid" scaffold in the list, run the matching protocol
	for ( iter_input_pdb_filename = pdb_file_names.begin(), iter_last_pdb = pdb_file_names.end(); iter_input_pdb_filename != iter_last_pdb; ++iter_input_pdb_filename ) {

		std::string scaffold_filename = *iter_input_pdb_filename;
		std::string scaffold_file_path = utility::file::FileName( *iter_input_pdb_filename ).path();
		std::string scaffold_file_basename = utility::file::FileName( *iter_input_pdb_filename ).base();

		if ( !utility::file::file_exists( scaffold_filename ) ) {
			std::cerr << "ERROR: Input pdb " << scaffold_filename << " not found." << std::endl;
			utility_exit();
		}

		// Read in the scaffold pdb with the aligned (rechained) epitope (i.e. the hybrid PDB file which has the epitope as chain Z and the rest of the scaffold)
		pose::Pose scaffold;
		import_pose::pose_from_pdb( scaffold, scaffold_filename );
		TR << "matching on scaffold '" << scaffold_filename << "'" << std::endl;

		// the option -scaffold_epitope_residue_range will have the stretch of residues (inclusive) that should be used to do the alignment. the remaining residues
		// in the scaffold should be scanned for positions which match the ones specified in -reference_structure_positions_to_match.

		if ( !option[ find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range ].user() ) {
			std::cerr << "ERROR: Epitope residue range on scaffold not specified." << std::endl;
			utility_exit();
		} else if ( option[ find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range ]().size() > 2 ) {
			std::cerr << "ERROR: Invalid epitope residue range on scaffold specified." << std::endl;
			utility_exit();
		}

		// need to convert the StringVectorOptionKey option to start and stop residues, because that's what we'll align on later.
		std::string scaffold_epitope_residue_start_pdb = option[ find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range ]()[ 1 ];
		std::string scaffold_epitope_residue_stop_pdb = option[ find_position_matches_using_inverse_rotamers::scaffold_epitope_residue_range ]()[ 2 ];

		TR << "scaffold epitope residue range (PDB numbering): " << scaffold_epitope_residue_start_pdb  << "-" << scaffold_epitope_residue_stop_pdb << std::endl;
			
		Size scaffold_epitope_residue_start;
		Size scaffold_epitope_residue_stop;

		if ( !option[ find_position_matches_using_inverse_rotamers::assume_command_line_has_pose_numbering ] ) {
			// convert PDB numbering to pose numbering, because we use the pose numbering later
			Size pdb_resnum; 
			char chain, icode;
			parse_position( scaffold_epitope_residue_start_pdb, chain, pdb_resnum, icode );
			scaffold_epitope_residue_start = (scaffold.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

			parse_position( scaffold_epitope_residue_stop_pdb, chain, pdb_resnum, icode );
			scaffold_epitope_residue_stop = (scaffold.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

		} else {
			std::istringstream iss( scaffold_epitope_residue_start_pdb );
			iss >> scaffold_epitope_residue_start;
			
			iss.clear();
			iss.str( scaffold_epitope_residue_stop_pdb );
			iss >> scaffold_epitope_residue_stop;
		}

		TR << "scaffold epitope residue range (pose numbering): " << scaffold_epitope_residue_start << "-" << scaffold_epitope_residue_stop << std::endl;

		// Make sure the number of epitope residues in the reference and scaffold structures is the same.
		Size number_epitope_residues_scaffold = scaffold_epitope_residue_stop - scaffold_epitope_residue_start + 1;
		Size number_epitope_residues_reference_structure = reference_structure_epitope_residue_stop - reference_structure_epitope_residue_start + 1;
		
		if ( number_epitope_residues_scaffold != number_epitope_residues_reference_structure ) {
			TR << "ERROR: reference epitope and scaffold epitope do not contain same number of residues!" << std::endl;
			utility_exit();
		}

		// align the scaffold to the reference structure using just the epitope residues in both.
		// use function superimpose_pose in rms_util.cc to do the alignment. 

		id::AtomID_Map< id::AtomID > atom_correspondence_map;
		pose::initialize_atomid_map( atom_correspondence_map, scaffold, id::BOGUS_ATOM_ID ); // maps every atomid to bogus atom
		
		// let's say the "match" on the scaffold is over residues 125-129, and the epitope in the reference structure is residues
		// 521-525.  if we iterate from 125 to 129, to get the offset we should use ii - scaffold_epitope_residue_start.
		Size reference_structure_residue;
		for ( Size ii = scaffold_epitope_residue_start; ii <= scaffold_epitope_residue_stop; ++ii ) {
			reference_structure_residue = reference_structure_epitope_residue_start + ( ii - scaffold_epitope_residue_start );
			id::AtomID const id1( scaffold.residue( ii ).atom_index("CA"), ii );
			id::AtomID const id2( reference_structure.residue( reference_structure_residue ).atom_index("CA"), reference_structure_residue );
			TR << "mapping reference structure residue " << reference_structure_residue << " to scaffold residue " << ii << std::endl;
			atom_correspondence_map[ id1 ] = id2;
		}
		Real alignment_rms = scoring::superimpose_pose( scaffold, reference_structure, atom_correspondence_map );
		TR << "scaffold aligned to reference_structure with rms '" << alignment_rms << "'" << std::endl;

		// For debugging, write out aligned structure
		std::string aligned_filename;
		if ( option[ find_position_matches_using_inverse_rotamers::output_aligned_scaffold ] ) {
			aligned_filename = "scaffold_aligned." +  scaffold_file_basename + ".pdb"; // already includes the ".pdb"
			TR << "Writing out aligned scaffold to file " << aligned_filename << std::endl;
			scaffold.dump_pdb( aligned_filename );
		}
		aligned_filename = scaffold_file_path + aligned_filename;

		// now we get to the part where we have to build inverse rotamers for each of the positions specified in the option
		// reference_structure_positions_to_match

		for ( Size ii = 1; ii <= positions_to_match.size(); ++ii ) {

			//Size const position = option[ find_position_matches_using_inverse_rotamers::reference_structure_positions_to_match ]()[ ii ];
			Size const position = positions_to_match[ ii ];
			conformation::Residue const & reference_res = reference_structure.residue( position );
			
			TR << "Building generic residues for residue " << reference_res.name3() << "-" << reference_res.seqpos() << std::endl;

			// use a function written for the enzyme design application to build residue specific rotamers
			core::chemical::ResidueTypeCOP invrot_restype( reference_res.type() );
			
			utility::vector1< conformation::ResidueOP > rotamers;
			build_bb_independent_rotamers( invrot_restype, rotamers );
			runtime_assert( rotamers.size() > 0 );
			TR << "built " << rotamers.size() << " bbindependent rotamers for residue " << rotamers[1]->type().name() << "." << std::endl;

			// so the function above should have built all possible rotamers of an asparagine residue. but these will all 
			// presumably be aligned on the backbone atoms (in an arbitrary coordinate system), with differing chi angles (obviously). 
			// what we need to do is align the ND2, OD1 and CG atoms of each rotamer to the reference structure residue so that the backbone
			// atoms of each rotamer are all over the place. once we have that we can go through every position in every scaffold and calculate
			// the rmsd (and other measurements) to every rotamer and see if there are any scaffolds that can accomodate the asn (and glycan) well.
			// not sure what the best function for aligning the inverse rotamer to the scaffold residue is. it seems like the enzyme design 
			// machinery is doing a manual alignment by using the set_xyz() function in Residue on every atom of the rotamer. it's precalculating
			// what the xyz values of every atom should be and then setting the inverse rotamer atom coordinates to those values. 

			// Invert the rotamer library as specified by the motif
			for( Size ir = 1; ir <= rotamers.size(); ++ir ) {
				
				// the Residue class contains alignment methods that will align one Residue onto another. specifically, there is a function which
				// takes an atom pair vector which is the set of atoms to use to do the alignment.
				// the function assumes 1) that the atom pairs are made up of strings of the atom type names, and 2) that the first atom pair
				// is of the "central" atom (though, I'm not sure how critical this is).
				utility::vector1< std::pair< std::string, std::string > > atom_pairs;
				atom_pairs.push_back( std::make_pair< std::string, std::string >( "CG", "CG" ) );
				atom_pairs.push_back( std::make_pair< std::string, std::string >( "ND2", "ND2" ) );
				atom_pairs.push_back( std::make_pair< std::string, std::string >( "OD1", "OD1" ) );
				
				rotamers[ ir ]->orient_onto_residue( reference_res, atom_pairs );

			}

#ifdef FILE_DEBUG
			TR << "Writing out inverse rotamers for residue " << reference_res.name3() << "-" << reference_res.seqpos() << std::endl;
			std::stringstream out;
			out << "inverse_rotamers." << scaffold_file_basename << "." << reference_res.name3() << "-" << reference_res.seqpos() << ".pdb";
			dump_rotamerset_pdb( rotamers, out.str() );
#endif

			// now, iterate over each position in the scaffold to see if any of the inverse rotamers are close to it.

			for ( Size ir = 1; ir <= rotamers.size(); ++ir ) {
			
				conformation::Residue const & inverse_rotamer = *(rotamers[ ir ]);

				Real chi1 = inverse_rotamer.chi( 1 ); if ( chi1 < 0 ) { chi1 += 360; } // Rosetta uses the bounds [-180,180] for all dihedral. convert to [0,360] range.
				Real chi2 = inverse_rotamer.chi( 2 ); if ( chi2 < 0 ) { chi2 += 360; }
				
				Real chi3 = 0.0; Real chi4 = 0.0;
				if ( option[ find_position_matches_using_inverse_rotamers::check_for_valid_chi3_chi4_dihedrals ] ) {
					chi3 = inverse_rotamer.chi( 3 ); if ( chi3 < 0 ) { chi3 += 360; }
					chi4 = inverse_rotamer.chi( 3 ); if ( chi4 < 0 ) { chi4 += 360; }
				}

//#ifdef FILE_DEBUG
				TR << "Calculated dihedral angles for inverse rotamer id " << ir << ", chi1: " << chi1 << ", chi2: " << chi2 << std::endl;
//#endif

				for ( Size scaffold_resnum = 1; scaffold_resnum <= scaffold.total_residue(); ++scaffold_resnum ) {

					conformation::Residue const & scaffold_res = scaffold.residue( scaffold_resnum );

#ifdef FILE_DEBUG
					TR << "Testing scaffold position " << scaffold_res.name3() << "-" << scaffold_res.seqpos() << std::endl;
#endif

					if ( !scaffold_res.has("CA") ) {
						TR << "Warning: position " << scaffold_res.seqpos() << " is missing CA atom. skipping this position." << std::endl;
						continue;
					}

					// skip residues that don't have a cbeta atom, unless it's a glycine. for glycines, use the 2HA atom position
					if ( scaffold_res.name3() != "GLY" && !scaffold_res.has( "CB" ) ) {
						TR << "Warning: position " << scaffold_res.seqpos() << " is missing CB atom, and is not a glycine. skipping this position." << std::endl;
						continue;
					}

					// use the rms_util function "rms_at_all_corresponding_atoms" to get the rmsd between the ca and cb atoms of the scaffold and rotamer.
					// actually, can't do that. all of the rms_util functions take a Pose object and the inverse rotamers are loose Residue objects.
					// since we're only doing ca and cb atoms, just do the calculation manually here.
					Real sum = 0.0;
					core::Vector ca_diff = inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ) - scaffold_res.xyz( scaffold_res.atom_index( "CA" ) );
					sum += ca_diff.length_squared();
					
					core::Vector cb_diff;
					if ( scaffold_res.name3() == "GLY" ) {
						// use the 2HA atom instead of CB for glycine
						cb_diff = inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ) - scaffold_res.xyz( scaffold_res.atom_index( "2HA" ) );

						// the best solution would be to determine what the coordinates of a cb atom would be if this residue were an alanine, using st like that below.
						//if ( pos.atomExists("N") && pos.atomExists("C") ) {
						//	CartesianPoint CBcoor = CartesianGeometry::build( pos.getAtom("CA").getCoor(), pos.getAtom("N").getCoor(), pos.getAtom("C").getCoor(), 1.521, 110.5, -122.5 );
						//	Atom cb( pos.getCurrentIdentity().getIdentityId() + ",CB", CBcoor );

					} else {
						cb_diff = inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ) - scaffold_res.xyz( scaffold_res.atom_index( "CB" ) );
					}
					sum += cb_diff.length_squared();
					
					Real rmsd = std::sqrt( sum / 2 );
					if ( rmsd < option[ find_position_matches_using_inverse_rotamers::ir_rmsd ] ) {

						if ( scaffold_res.name3() != "GLY" ) {
							TR << "calculated rmsd for point irCA: (" << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).x() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).y() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).z() 
								<< "), irCB: (" << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).x() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).y() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).z() 
								<< "), sCA: (" << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).x() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).y() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).z()
								<< "), sCB: (" << scaffold_res.xyz( scaffold_res.atom_index( "CB" ) ).x() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CB" ) ).y() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CB" ) ).z()
								<< ")" << std::endl;
						} else {
							TR << "calculated rmsd for point irCA: (" << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).x() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).y() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CA" ) ).z() 
								<< "), irCB: (" << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).x() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).y() << ", " << inverse_rotamer.xyz( inverse_rotamer.atom_index( "CB" ) ).z() 
								<< "), sCA: (" << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).x() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).y() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "CA" ) ).z()
								<< "), sCB: (" << scaffold_res.xyz( scaffold_res.atom_index( "2HA" ) ).x() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "2HA" ) ).y() << ", " << scaffold_res.xyz( scaffold_res.atom_index( "2HA" ) ).z()
								<< ")" << std::endl;
						}

						TR << "scaffold position " << scaffold_resnum << " has rmsd: " << ObjexxFCL::format::F( 4,2,rmsd ) << " to inverse rotamer: " << ir << ". checking for glycan clashes." << std::endl;

						// also scan the hit for glycan clashes
						Size number_glycan_clashes = 0;
						
						for ( Size glycan_res_index = 1; glycan_res_index <= reference_structure_glycan_residues.size(); ++glycan_res_index ) {
							conformation::Residue const & glycan_res = reference_structure.residue( reference_structure_glycan_residues[ glycan_res_index ] );

							for ( Size glycan_atom_index = 1; glycan_atom_index <= glycan_res.natoms(); ++glycan_atom_index ) {
								
								for ( Size scaffold_res_index = 1; scaffold_res_index <= scaffold.total_residue(); ++scaffold_res_index ) {
									conformation::Residue const & scaffold_res = scaffold.residue( scaffold_res_index );

									for ( Size scaffold_atom_index = 1; scaffold_atom_index <= scaffold.residue( scaffold_res_index ).natoms(); ++scaffold_atom_index ) {
									
										Real dist_sq = scaffold_res.xyz( scaffold_atom_index ).distance_squared( glycan_res.xyz( glycan_atom_index ) );
										if ( dist_sq < 4.0 ) {
#ifdef FILE_DEBUG
											TR << "scaffold residue '" << scaffold_res_index << "' atom '" << scaffold_res.atom_name( scaffold_atom_index ) 
												<< "' comes within 4.0 Ang of glycan residue '" << reference_structure_glycan_residues[ glycan_res_index ]
												<< "' atom '" << glycan_res.atom_name( scaffold_atom_index ) << "'. dist_sq: " << ObjexxFCL::format::F( 4,2,dist_sq ) << std::endl;
#endif
											number_glycan_clashes++;
										}

									} // end scaffold residue atoms

								} // end scaffold residues
							
							} // end glycan residue atoms
							
						} // end glycan residues

						TR << "\tnumber of glycan clashes " << number_glycan_clashes << " beneath maximum of " << number_glycan_clashes_allowed << ". checking dihedral angles." << std::endl;

						if ( number_glycan_clashes <= number_glycan_clashes_allowed ) {

							// do one more filtering step: check that the dihedral angles are in regions observed for glycosylated ASNs (as described in Petrescu, Glycobiology, 2004)
							// chi1 must be either between [260,320] or [165,225]
							if ( !( ( 260 < chi1 && chi1 < 320 ) || ( 165 < chi1 && chi1 < 225 ) ) ) {
								TR << "\t\tfailed chi1 check. chi1 must be either between [260,320] or [165,225]. chi1 value was " << chi1 << std::endl;
								//continue;
							}
							// chi2 must be between [90,260]
							if ( !( 90 < chi2 && chi2 < 260 ) ) {
								TR << "\t\tfailed chi2 check. chi2 must be between [90,260]. chi2 value was " << chi2 << std::endl;
								//continue;
							}

							if ( option[ find_position_matches_using_inverse_rotamers::check_for_valid_chi3_chi4_dihedrals ] ) {
								// chi3 must be between 176.8 +/- 12.0, so [164,188]
								if ( !( 164 < chi3 && chi3 < 188 ) ) { continue; }
								// chi4 must be between 260.2 +/- 21.1, so [239.1,281.3]
								if ( !( 239.1 < chi4 && chi4 < 281.3 ) ) { continue; }
							}
							
							// passed all dihedral checks...
							/*VectorPair bond_vectors( scaffoldCaCb(0).getCoor(), scaffoldCaCb(1).getCoor(), inverseCaCb(0).getCoor(), inverseCaCb(1).getCoor() );
							//TR << "a1/scCA: " << scaffoldCaCb(0).getCoor()
							//				<< ", a2/scCB: " << scaffoldCaCb(1).getCoor()
							//				<< ", b1/irCA: " << inverseCaCb(0).getCoor()
							//				<< ", b2/irCB: " << inverseCaCb(1).getCoor() << std::endl;
							//bond_vectors.calcAll();
							bond_vectors.calcDotProductAngle();
							double dotp_angle = bond_vectors.getDotProductAngle();
							bond_vectors.calcTorsion();
							double dihedral = bond_vectors.getTorsion();
							*/
							
							if ( option[ find_position_matches_using_inverse_rotamers::check_for_valid_chi3_chi4_dihedrals ] ) {
								//fprintf( fout, "%c/%d%c\t%d\t%s\t%c/%d%c\t%4.2f\t%d\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%4.1f\n",
								fprintf( fout, "%c/%d%c\t%d\t%s\t%s\t%d\t%c/%c%d%c\t%5.3f\t%d\t%4.1f\t%4.1f\t%4.1f\t%4.1f\n",
									reference_structure.pdb_info()->chain( position ),
									(int)reference_structure.pdb_info()->number( position ),
									reference_structure.pdb_info()->icode( position ),
									(int)ir,
									scaffold_filename.c_str(),
									aligned_filename.c_str(),
									(int)scaffold.total_residue(),
									scaffold.pdb_info()->chain( scaffold_resnum ),
									scaffold.residue( scaffold_resnum ).name1(),
									(int)(scaffold.pdb_info()->number( scaffold_resnum )),
									scaffold.pdb_info()->icode( scaffold_resnum ),
									rmsd,
									(int)number_glycan_clashes,
									chi1, chi2, chi3, chi4 );
									//, dotp_angle, dihedral );
							} else {
								//fprintf( fout, "%c/%d%c\t%d\t%s\t%c/%d%c\t%4.2f\t%d\t%4.1f\t%4.1f\t%4.1f\t%4.1f\n",
								fprintf( fout, "%c/%d%c\t%d\t%s\t%s\t%d\t%c/%c%d%c\t%5.3f\t%d\t%4.1f\t%4.1f\n",
									reference_structure.pdb_info()->chain( position ),
									(int)reference_structure.pdb_info()->number( position ),
									reference_structure.pdb_info()->icode( position ),
									(int)ir,
									scaffold_filename.c_str(),
									aligned_filename.c_str(),
									(int)scaffold.total_residue(),
									scaffold.pdb_info()->chain( scaffold_resnum ),
									scaffold.residue( scaffold_resnum ).name1(),
									(int)(scaffold.pdb_info()->number( scaffold_resnum )),
									scaffold.pdb_info()->icode( scaffold_resnum ),
									rmsd,
									(int)number_glycan_clashes,
									chi1, chi2 );
									//, dotp_angle, dihedral );
							}
							fflush( fout );

							TR << "\tall checks passed. writing out match." << std::endl;

							// write out inverse rotamers to check that they were done properly
							if ( option[ find_position_matches_using_inverse_rotamers::output_inverse_rotamers ] ) {
								std::stringstream ss;
								ss << "ir." << scaffold_file_basename << "." << position << "-" << ir << ".pdb";
								std::ofstream out( ss.str().c_str() );
								Size atom_counter = 0;
								core::io::pdb::dump_pdb_residue( inverse_rotamer, atom_counter, out );
								out.close();
							}	
						
						} else {
							TR << "\t\tscaffold postion " << scaffold_resnum << " failed clash check. number_glycan_clashes: " << number_glycan_clashes
								<< ", number_glycan_clashes_allowed: " << number_glycan_clashes_allowed << std::endl;
						}

#ifdef FILE_DEBUG
					} else {
						TR << "scaffold position " << scaffold_resnum << " inverse rotamer: " << ir << " failed rmsd check. ca_diff.length_squared: "
							<< ca_diff.length_squared() << ", cb_diff.length_squared(): " << cb_diff.length_squared() 
							<< ", rmsd: " << ObjexxFCL::format::F( 4,2,rmsd )
							<< ", cutoff: " << option[ find_position_matches_using_inverse_rotamers::ir_rmsd ] << std::endl;
#endif
					}

				} // end loop over scaffold positions

			} // end loop over inverse rotamers

		} // end loop over positions to look for matches
	}

	// close the results output stream
	fclose( fout );
	TR << "Done." << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

