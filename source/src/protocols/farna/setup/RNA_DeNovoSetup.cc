// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/setup/RNA_DeNovoSetup.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/setup/RNA_DeNovoSetup.hh>
#include <protocols/farna/setup/RNA_DeNovoParameters.hh>
#include <protocols/farna/secstruct/RNA_SecStruct.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/util.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/leontis_westhof_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>

static basic::Tracer TR( "protocols.farna.setup.RNA_DeNovoSetup" );
using utility::vector1;
using utility::tools::make_vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::A;
using namespace core;
using namespace core::chemical::rna;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
// @details
//
// Setup of rna_denovo app, or (fragment assembly of RNA, FARNA) protocol
//
// Takes command-line input ( -fasta, -secstruct_file, -obligate_pair, etc.) and
//  produces objects needed for RNA_DeNovoProtocol:
//
//    pose  (starter extended pose, with constraints, rna_data, full_model_info, etc.)
//    native_pose
//    rna_de_novo_parameters (stems in secondary structure, obligate pairs, etc.)
//    refine_pose_list (if in -refine_silent_file mode).
//
// Refactored setup stuff out of rna_denovo.cc
// Ported over tools/rna_tools/bin/rna_denovo_setup.py
//       -- Rhiju Das, June 2016
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace farna {
namespace setup {

//Constructor
RNA_DeNovoSetup::RNA_DeNovoSetup()
{}

//Destructor
RNA_DeNovoSetup::~RNA_DeNovoSetup()
{}

void
RNA_DeNovoSetup::initialize_from_command_line()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace protocols::farna::options;

	options_ = RNA_DeNovoProtocolOptionsOP( new RNA_DeNovoProtocolOptions);
	options_->initialize_from_command_line();

	rsd_set_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	if ( options_->use_legacy_setup() || options_->rna_params_file().size() > 0 ) {
		de_novo_setup_from_command_line_legacy();
	} else {
		de_novo_setup_from_command_line();
	}

	// if output_res_num supplied, this will change PDBInfo numbering & chain.
	set_output_res_and_chain( *pose_, option[ OptionKeys::rna::farna::output_res_num ].resnum_and_chain() );

	// refine_pose is a seldom-used functionality at the moment -- not well tested.
	setup_refine_pose_list();

}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
// Following is directly adapted from rna_denovo_setup.py
//
// The main update is that for most of the script, we previously
//  kept track of 'conventional' numbers & chains for most
//  residues (working_res, cutpoint_open, extra_minimize_res, etc.)
//
// Now most of those values are held in 'full model' numbering.
//
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
void
RNA_DeNovoSetup::de_novo_setup_from_command_line()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::farna::options;
	using namespace protocols::farna::secstruct;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	////////////////////
	// Step 1
	////////////////////
	vector1< std::string > sequence_strings  = option[ OptionKeys::rna::farna::sequence ]();
	vector1< std::string > fasta_files = option[ in::file::fasta ]();
	int const offset = option[ OptionKeys::rna::farna::offset ]();

	// Sequence setup:
	// FullModelParameters is a nice object that holds sequence, non-standard residues ("Z[Mg]"),
	// and what chains and residue numbers to use.
	FullModelParametersOP full_model_parameters;
	vector1< Size > cutpoint_open_in_full_model;
	if ( fasta_files.size() > 0 ) {
		// use fasta readin developed for stepwise application -- also reads in
		// numbers & chains based on fasta header lines.
		runtime_assert( sequence_strings.size() == 0 );
		runtime_assert( fasta_files.size() == 1 );
		full_model_parameters = stepwise::setup::get_sequence_information( fasta_files[ 1 ], cutpoint_open_in_full_model );
		if ( offset != 0 ) {
			vector1< int > new_numbering = full_model_parameters->conventional_numbering();
			for ( Size n = 1; n <= new_numbering.size(); n++ ) { new_numbering[ n ] += offset; }
			full_model_parameters->set_conventional_numbering( new_numbering );
		}
	} else {
		// basic read-in of sequence from command line
		runtime_assert( sequence_strings.size() > 0 );
		std::string sequence( sequence_strings[1] );
		for ( Size n = 2; n <= sequence_strings.size(); n++ ) sequence += std::string( " " + sequence_strings[ n ] );
		cutpoint_open_in_full_model = core::sequence::strip_spacers( sequence );
		std::map< Size, std::string > non_standard_residue_map = stepwise::setup::parse_out_non_standard_residues( sequence );
		vector1< int > res_numbers_in_pose;
		for ( Size n = 1; n <= sequence.size(); n++ ) res_numbers_in_pose.push_back( n + offset );
		stepwise::setup::get_extra_cutpoints_from_names( sequence.size(), cutpoint_open_in_full_model, non_standard_residue_map );
		full_model_parameters = FullModelParametersOP( new FullModelParameters( sequence, cutpoint_open_in_full_model, res_numbers_in_pose ) );
		full_model_parameters->set_non_standard_residue_map( non_standard_residue_map );
	}
	std::string const sequence = full_model_parameters->full_sequence();

	////////////////////
	// Step 2
	////////////////////
	// Other useful residues.
	if ( option[ full_model::cutpoint_open ].user() ) {
		cutpoint_open_in_full_model  =
			full_model_parameters->conventional_to_full( option[ full_model::cutpoint_open ].resnum_and_chain() );
	}
	vector1< Size > working_res        =
		full_model_parameters->conventional_to_full( option[ full_model::working_res ].resnum_and_chain() ); //all working stuff
	vector1< Size > const cutpoint_closed          =
		full_model_parameters->conventional_to_full( option[ full_model::cutpoint_closed ].resnum_and_chain() );
	vector1< Size > extra_minimize_res =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::minimize::extra_minimize_res ].resnum_and_chain() );
	vector1< Size > extra_minimize_chi_res =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::minimize::extra_minimize_chi_res ].resnum_and_chain() );
	vector1< Size > input_res_user_defined =
		full_model_parameters->conventional_to_full( option[ in::file::input_res ].resnum_and_chain() );
	vector1< Size > input_silent_res_user_defined =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::input_silent_res ].resnum_and_chain() );
	vector1< Size > virtual_anchor =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::virtual_anchor ].resnum_and_chain() );
	vector1< Size > obligate_pair =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::obligate_pair ].resnum_and_chain() );
	vector1< Size > remove_pair =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::remove_pair ].resnum_and_chain() );
	vector1< Size > remove_obligate_pair =
		full_model_parameters->conventional_to_full( option[ OptionKeys::rna::farna::remove_obligate_pair ].resnum_and_chain() );

	////////////////////
	// Step 3
	////////////////////
	// secondary structure setup.
	RNA_SecStruct secstruct( option[ OptionKeys::rna::farna::secstruct ](), option[ OptionKeys::rna::farna::secstruct_file ](), sequence );
	// "general" secondary structure includes non-canonical pairs that should be connected by jumps during run; used with -bps_moves.
	RNA_SecStruct secstruct_general( option[ OptionKeys::rna::farna::secstruct_general ](), option[ OptionKeys::rna::farna::secstruct_general_file ](), sequence );

	secstruct.check_compatible_with_sequence( sequence, true  /*check_complementarity*/ );
	secstruct_general.check_compatible_with_sequence( sequence, false /*check_complementarity*/ );

	if ( !secstruct_general.blank() && !options_->bps_moves() ) utility_exit_with_message("cannot supply secstruct_general without bps_moves");

	////////////////////
	// Step 4
	////////////////////
	/////////////
	vector1< std::string > input_pdbs = option[ in::file::s ]();
	vector1< std::string > input_silent_files = option[ in::file::silent ]();
	std::string const working_native_pdb = option[ OptionKeys::rna::farna::working_native ]();
	vector1< std::string > obligate_pair_explicit = option[ OptionKeys::rna::farna::obligate_pair_explicit ]();
	vector1< std::string > chain_connections = option[ OptionKeys::rna::farna::chain_connection ]();

	runtime_assert( remove_pair.size() % 2 == 0 );
	for ( Size n = 1; n <= remove_pair.size(); n += 2 ) {
		secstruct.remove_pair( std::make_pair( remove_pair[ n ], remove_pair[ n + 1 ] ) );
	}

	////////////////////
	// Step 5
	////////////////////
	// in full_model numbering.
	vector1< vector1< int > > resnum_list;
	vector1< Size > input_res;
	Size input_res_user_defined_count( 0 );
	for ( Size n = 1; n <= input_pdbs.size(); n++ ) {
		std::string const pdb = input_pdbs[ n ];
		std::string pdb_seq;
		vector1< int > resnum;
		vector1< char >  chain;
		get_seq_and_resnum( pdb, pdb_seq, resnum, chain);
		vector1< Size > resnum_in_full_model;

		if ( input_res_user_defined_count + resnum.size() <= input_res_user_defined.size() ) {
			// input res could have come from user after flag -input_res
			for ( Size q = 1; q <= resnum.size(); q++ ) {
				input_res_user_defined_count++;
				resnum_in_full_model.push_back( input_res_user_defined[ input_res_user_defined_count ] );
			}
		} else {
			// figure out residue numbers from PDB resnum & chain.
			resnum_in_full_model = full_model_parameters->conventional_to_full( std::make_pair( resnum, chain ) );
		}

		std::string actual_seq = "";
		for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << resnum[q] << " " << chain[q] << " exists in two pdb files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res.push_back( i );
		}
		if ( pdb_seq != actual_seq ) {
			TR << TR.Red << pdb_seq << std::endl;
			TR << TR.Red << actual_seq << std::endl;
			utility_exit_with_message("The sequence in "+pdb+" does not match input sequence!!");
		}
		resnum_list.push_back( resnum_in_full_model );
	}


	////////////////////
	// Step 6
	////////////////////
	// go through silent files for each residue, and sequences match up.
	Size input_silent_res_user_defined_count( 0 );
	vector1< Size > input_silent_res;
	if ( input_silent_res_user_defined.size() > 0 ) {
		// must have run through input_res.
		runtime_assert( input_res_user_defined_count == input_res_user_defined.size() );
		input_silent_res = input_silent_res_user_defined;
	} else {
		input_silent_res = input_res_user_defined;
		input_silent_res_user_defined_count = input_res_user_defined_count;
	}

	for ( Size n = 1; n <= input_silent_files.size(); n++ ) {
		std::string const silent = input_silent_files[ n ];
		std::string seq = get_silent_seq( silent );
		Size len_seq = seq.size();
		vector1< Size > resnum_in_full_model;
		if ( ( input_silent_res_user_defined_count + len_seq ) <= input_silent_res.size() ) {
			// input res could have come from user after flag -input_res or -input_silent_res
			for ( Size q = 1; q <= len_seq; q++ ) {
				input_silent_res_user_defined_count++;
				resnum_in_full_model.push_back( input_silent_res[ input_silent_res_user_defined_count ] );
			}
		} else {
			///////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////
			// would need to read in working sequence used in silent file.
			///////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////
		}
		runtime_assert( resnum_in_full_model.size() == len_seq);

		std::string actual_seq;
		for ( Size q = 1; q <= len_seq; q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << i << " exists in two pdb/silent files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res.push_back( i );
		}

		if ( seq != actual_seq ) {
			TR << TR.Red << seq << std::endl;
			TR << TR.Red << actual_seq << std::endl;
			utility_exit_with_message("The sequence in "+silent+" does not match input sequence!!");
		}

		resnum_list.push_back( resnum_in_full_model );
	}
	runtime_assert( input_silent_res_user_defined_count == input_silent_res.size() );


	////////////////////
	// Step 7
	////////////////////
	runtime_assert( obligate_pair_explicit.size() % 5 == 0 );
	vector1< Size > obligate_pair_explicit_full_model;
	for ( Size m = 0; m < obligate_pair_explicit.size()/5; m++ ) {
		std::vector< int > resnum;
		std::vector< char > chains;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains );
		resnum_full = full_model_parameters->conventional_to_full( std::make_pair( vector1<int>( resnum ),
			vector1<char>( chains ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains );
		resnum_full = full_model_parameters->conventional_to_full( std::make_pair( vector1<int>( resnum ),
			vector1<char>( chains ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos2 = resnum_full[ 1 ];

		obligate_pair_explicit_full_model.push_back( pos1 );
		obligate_pair_explicit_full_model.push_back( pos2 );
	}


	vector1< char > const & conventional_chains = full_model_parameters->conventional_chains();
	vector1< Size > domain_map( sequence.size(), 0 );
	// Go through each of the inputs, and look for broken chains -- will define obligate pairs.
	for ( Size n = 1; n <= resnum_list.size(); n++ ) {
		vector1< Size > resnum = resnum_list[ n ];
		//Find obligate pairs
		vector1< vector1< Size > > chunks;
		vector1< Size > curr_chunk;
		Size i( 0 ), j( 0 );
		char c(' '), d( ' ' );
		for ( Size q = 1; q <= resnum.size(); q++ ) {
			i = resnum[ q ];
			c = conventional_chains[ resnum[ q ] ];
			domain_map[ i ] = n;
			if ( j > 1 &&
					( ( ( i - 1 ) != j ) || d != c || cutpoint_open_in_full_model.has_value( j )  ) ) {
				chunks.push_back(curr_chunk);
				curr_chunk.clear();
			}
			j = i;
			d = c;
			curr_chunk.push_back( i );
		}
		chunks.push_back(curr_chunk);

		Size n_jumps( 0 );
		// look to see if each chunk might already be connected to the next.
		vector1< bool > connected_to_next( chunks.size(), false );
		for ( Size i = 1; i <= chunks.size(); i++ ) {
			Size const i_next = ( i < chunks.size() ) ? ( i + 1 ) : 1;

			// first check if there is already an obligate_pair between these chunks:
			bool found_pair( false );
			for ( Size j = 1; j <= chunks[ i ].size(); j++ ) {
				for ( Size k = 1; k <= chunks[ i_next  ].size(); k++ ) {
					Size const segment1_end   = chunks[ i     ][ j ];
					Size const segment2_start = chunks[ i_next ][ k ];
					Size const new_pos1 = std::min( segment1_end, segment2_start );
					Size const new_pos2 = std::max( segment1_end, segment2_start );
					vector1< Size > new_pair = make_vector1( new_pos1, new_pos2 );
					if ( already_listed_in_obligate_pair( new_pair, obligate_pair, obligate_pair_explicit_full_model ) ) {
						found_pair = true; break;
					}
					if ( found_pair ) break;
				}
			}
			if ( found_pair ) {
				connected_to_next[ i ] = true;
				n_jumps += 1;
			}
		}

		for ( Size i = 1; i <= chunks.size(); i++ ) {
			//obligate pairs -- go from end of one chain to beginning of next chain.
			if ( n_jumps == chunks.size() - 1 ) break;
			if ( connected_to_next[ i ] ) continue;
			Size const i_next = ( i < chunks.size() ) ? ( i + 1 ) : 1;

			Size const segment1_end   = chunks[ i     ][ chunks[ i ].size() ];
			Size const segment2_start = chunks[ i_next ][ 1 ];
			Size const new_pos1 = std::min( segment1_end, segment2_start );
			Size const new_pos2 = std::max( segment1_end, segment2_start );
			obligate_pair.push_back( new_pos1 );
			obligate_pair.push_back( new_pos2 );
			//   TR << "Creating new obligate pair: " << obligate_pair << " for chunk with residues " << resnum << std::endl;
			n_jumps++;
		}
	}

	vector1< std::pair< Size, Size > > canonical_pairs = flatten( secstruct.get_all_stems() );
	vector1< std::pair< Size, Size > > general_pairs   = flatten( secstruct_general.get_all_stems() );
	for ( Size n = 1; n <= general_pairs.size(); n++ ) {
		std::pair< Size, Size > const & p = general_pairs[ n ];
		if ( !canonical_pairs.has_value( p ) &&
				( domain_map[ p.first ] == 0 || domain_map[ p.second ] == 0 ) ) {
			vector1< Size > const new_pair = make_vector1( p.first, p.second );
			if ( !already_listed_in_obligate_pair( new_pair, obligate_pair, obligate_pair_explicit_full_model ) ) {
				obligate_pair.push_back( p.first );
				obligate_pair.push_back( p.second );
			}
		}
	}

	////////////////////
	// Step 8
	////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out working sequence, secstruct, secstruct_general
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( working_res.size() == 0 ) {
		for ( Size n = 1; n <= sequence.size(); n++ ) working_res.push_back( n );
	}

	// previously, kept spaces in sequences to mimic user-input with spaces.
	// don't worry about that anymore -- try to handle above in get_sequence_info.
	std::string working_sequence = working_res_map( sequence, working_res );
	RNA_SecStruct working_secstruct         = working_res_map( secstruct, working_res );
	RNA_SecStruct working_secstruct_general = working_res_map( secstruct_general, working_res );

	TR << "Sequence:            " << working_sequence << std::endl;
	TR << "Secstruct:           " << working_secstruct.secstruct() << std::endl;
	TR << "Secstruct [general]: " << working_secstruct_general.secstruct() << std::endl;

	////////////////////
	// Step 9
	////////////////////
	///////////////////////////////////////////////////////////////////
	// initialize variables needed for RNA_DeNovoParams (rna_params_)
	///////////////////////////////////////////////////////////////////
	vector1< Size > working_cutpoint_open   = working_res_map( cutpoint_open_in_full_model, working_res );
	vector1< Size > working_cutpoint_closed = working_res_map( cutpoint_closed, working_res );
	vector1< Size > working_virtual_anchor  = working_res_map( virtual_anchor, working_res );
	vector1< pose::rna::BasePair > working_obligate_pairs;
	vector1< vector1< std::pair< Size, Size > > > working_stems;

	for ( Size n = 1; n < working_res.size(); n++ ) {
		if ( working_res[ n+1 ] > working_res[n] + 1  && !working_cutpoint_open.has_value( n ) )  working_cutpoint_open.push_back( n );
	}

	/////////////////////////
	// working stems
	/////////////////////////
	working_stems = working_secstruct.get_all_stems( working_sequence, working_cutpoint_open );
	vector1< Size > working_input_res = working_res_map( input_res, working_res );
	if ( options_->fixed_stems() ) {
		update_working_obligate_pairs_with_stems( working_obligate_pairs, working_stems, working_input_res );
	}

	////////////////////
	// Step 10
	////////////////////
	////////////////////////
	// working pose
	////////////////////////
	Pose full_pose;
	pose_ = PoseOP( new Pose );
	std::string const full_annotated_sequence = full_model_parameters->full_annotated_sequence();
	make_pose_from_sequence( full_pose, full_annotated_sequence, *rsd_set_ );
	set_output_res_and_chain( full_pose, std::make_pair( full_model_parameters->conventional_numbering(),
		full_model_parameters->conventional_chains() ) );
	pdbslice( *pose_, full_pose, working_res );

	////////////////////
	// Step 11
	////////////////////
	////////////////////////
	// working constraints.
	////////////////////////
	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set( new ConstraintSet );
		option[ OptionKeys::constraints::force_pdb_info_mapping ].def( true ); // using option as global variable due to difficulty in dealing with static functions.
		cst_set = ConstraintIO::get_instance()->read_constraints( option[ OptionKeys::constraints::cst_file ](1), ConstraintSetOP( new ConstraintSet ), full_pose );
		full_pose.constraint_set( cst_set );
		id::SequenceMappingOP sequence_map( new id::SequenceMapping( working_res ) );
		sequence_map->reverse();
		pose_->constraint_set( full_pose.constraint_set()->remapped_clone( full_pose, *pose_, sequence_map ) );
	}

	if ( options_->cst_gap() ) {
		// also create a constraints file that will help close chains across 2-5 residue gaps...
		// this is a slight hack, but if it works, might be worth putting a term into Rosetta, as
		// well as automated handling of "working_res"
		for ( Size m = 1; m <= working_res.size(); m++ ) {
			Size const gap = working_res[m+1] - working_res[m];
			if ( gap > 1 && gap < 6 ) {
				Distance const stdev( 10.0 );
				Distance max_dist = gap * 5.0 + 4;
				Real const bonus = 200.0;
				//        cst_file_outstring +=  " O3* %d  C5* %d   FADE %6.3f  %6.3f  %6.3f %6.3f %6.3f \n" %
				//            ( cst_gap[0], cst_gap[1],  -stdev, max_dist, stdev, -1*bonus, bonus)
				FuncOP distance_func( new FadeFunc( -stdev, max_dist, stdev, -1.0*bonus, bonus ) );
				pose_->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint(
					named_atom_id_to_atom_id( NamedAtomID( " O3'",m), *pose_),
					named_atom_id_to_atom_id( NamedAtomID( " C5'",m+1), *pose_),
					distance_func ) ) ) );
			}
		}
	}

	////////////////////
	// Step 12
	////////////////////
	////////////////////////
	// working data
	////////////////////////
	std::string const in_path = option[ in::path::path ]()[1];
	if ( option[ OptionKeys::rna::data_file].user() ) {
		core::io::rna::RNA_DataReader rna_data_reader( in_path + option[ OptionKeys::rna::data_file ]  );
		// note that this actually does look at conventional numbering in a smart way (but not chains yet):
		rna_data_reader.fill_rna_data_info( *pose_ );
	}

	////////////////////
	// Step 13
	////////////////////
	////////////////////////////
	// working_obligate_pairs.
	////////////////////////////
	runtime_assert( obligate_pair.size() % 2 == 0 );
	for ( Size m = 0; m < obligate_pair.size()/2; m++ ) {
		Size const pos1 = obligate_pair[ 2*m + 1 ];
		Size const pos2 = obligate_pair[ 2*m + 2 ];
		if ( !working_res.has_value( pos1 ) ) continue;
		if ( !working_res.has_value( pos2 ) ) continue;

		bool rm_pair( false );
		for ( Size q = 1; q <= remove_obligate_pair.size(); q++ ) {
			Size const rm_pos1 = remove_obligate_pair[ 2*q + 1 ];
			Size const rm_pos2 = remove_obligate_pair[ 2*q + 2 ];
			if ( rm_pos1 == pos1 && rm_pos2 == pos2 ) {
				rm_pair = true; break;
			}
			if ( rm_pos1 == pos2 && rm_pos2 == pos1 ) {
				rm_pair = true; break;
			}
		}
		if ( rm_pair ) continue;

		working_obligate_pairs.push_back( BasePair( working_res.index( pos1 ), working_res.index( pos2 ), ANY_BASE_EDGE, ANY_BASE_EDGE, ANY_BASE_DOUBLET_ORIENTATION ) );
	}

	////////////////////
	// Step 14
	////////////////////
	//////////////////////////////////////
	// working_obligate_pairs [explicit]
	//////////////////////////////////////
	runtime_assert( obligate_pair_explicit.size() % 5 == 0 );
	for ( Size m = 0; m < obligate_pair_explicit.size()/5; m++ ) {

		std::vector< int > resnum;
		std::vector< char > chains;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains );
		resnum_full = full_model_parameters->conventional_to_full( std::make_pair( vector1<int>( resnum ),
			vector1<char>( chains ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];
		if ( !working_res.has_value( pos1 ) ) continue;

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains );
		resnum_full = full_model_parameters->conventional_to_full( std::make_pair( vector1<int>( resnum ),
			vector1<char>( chains ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos2 = resnum_full[ 1 ];
		if ( !working_res.has_value( pos2 ) ) continue;

		runtime_assert( obligate_pair_explicit[ 5*m + 3 ].size() == 1 );
		BaseEdge const edge1 = get_edge_from_char( obligate_pair_explicit[ 5*m + 3 ][ 0 ] );

		runtime_assert( obligate_pair_explicit[ 5*m + 4 ].size() == 1 );
		BaseEdge const edge2 = get_edge_from_char( obligate_pair_explicit[ 5*m + 4 ][ 0 ] );

		runtime_assert( obligate_pair_explicit[ 5*m + 5 ].size() == 1 );
		char const o_char = obligate_pair_explicit[ 5*m + 5 ][ 0 ];

		BaseDoubletOrientation orientation;
		if ( o_char == 'C' || o_char == 'T' ) { // convert from leontis-westhof cis/trans to antiparallel/antiparallel
			orientation  = get_base_doublet_orientation_from_LW( edge1, edge2, get_LW_orientation_from_char( o_char ) );
		} else {
			orientation = get_orientation_from_char( o_char );
		}

		working_obligate_pairs.push_back( BasePair( working_res.index( pos1 ), working_res.index( pos2 ),
			edge1, edge2, orientation ) );
	}

	////////////////////
	// Step 15
	////////////////////
	//////////////////////////////
	// working chain connections
	//////////////////////////////
	vector1< std::pair< vector1< Size >, vector1< Size > > > working_chain_connections;
	if ( chain_connections.size() > 0 ) {
		runtime_assert( chain_connections.size() >= 4 );
		if ( chain_connections.has_value( "SET1" ) ) {
			// new format is more flexible enough to allow for noncontiguous numbering in set1 vs. set2)
			Size which_set = 0;
			vector1< Size > resnum1;
			vector1< Size > resnum2;
			for ( Size k = 1; k <= chain_connections.size() + 1; k++ ) {
				if ( k == ( chain_connections.size() + 1 ) || chain_connections[ k ] == "SET1" ) { // at of previous block
					vector1< Size > working_resnum1 = working_res_map( resnum1, working_res );
					vector1< Size > working_resnum2 = working_res_map( resnum2, working_res );
					if ( working_resnum1.size() > 0 &&
							working_resnum2.size() > 0 ) {
						working_chain_connections.push_back( std::make_pair( working_resnum1, working_resnum2) );
					}
					resnum1.clear();
					resnum2.clear();
					which_set = 1;
					continue;
				}
				if ( chain_connections[ k ] == "SET2" ) {
					which_set = 2;
					continue;
				}
				runtime_assert( which_set > 0 );

				std::vector< int > resnum;
				std::vector< char > chains;
				bool ok = utility::get_resnum_and_chain_from_one_tag( chain_connections[ k ], resnum, chains );
				runtime_assert( ok );
				vector1< Size > resnum_full = full_model_parameters->conventional_to_full(
					std::make_pair( vector1<Size>( resnum ),
					vector1< char >( chains ) ) );
				vector1< Size > & resnum_for_set = ( which_set == 1 ) ? resnum1 : resnum2;
				for ( Size n = 1; n <= resnum_full.size(); n++ ) resnum_for_set.push_back( resnum_full[ n ] );
			}
		} else { // legacy format
			runtime_assert( chain_connections.size() % 4 == 0 );
			Size const n_connect = chain_connections.size() / 4;
			for ( Size i = 0; i < n_connect; i++ ) {
				Size const curr_0 = i * 4;
				Size const seg1_start = ObjexxFCL::int_of( chain_connections[curr_0 + 1] );
				Size const seg1_stop  = ObjexxFCL::int_of( chain_connections[curr_0 + 2] );
				Size const seg2_start = ObjexxFCL::int_of( chain_connections[curr_0 + 3] );
				Size const seg2_stop  = ObjexxFCL::int_of( chain_connections[curr_0 + 4] );
				if ( !working_res.has_value( seg1_start ) ) continue;
				if ( !working_res.has_value( seg1_stop  ) ) continue;
				if ( !working_res.has_value( seg2_start ) ) continue;
				if ( !working_res.has_value( seg2_stop ) ) continue;

				vector1< Size > working_resnum1, working_resnum2;
				for ( Size n = seg1_start; n <= seg1_stop; n++ ) {
					runtime_assert( working_res.has_value( n ) );
					working_resnum1.push_back( working_res.index( n ) );
				}
				for ( Size n = seg2_start; n <= seg2_stop; n++ ) {
					runtime_assert( working_res.has_value( n ) );
					working_resnum2.push_back( working_res.index( n ) );
				}
				working_chain_connections.push_back( std::make_pair( working_resnum1, working_resnum2) );
			}
		}
	}

	////////////////////
	// Step 16
	////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	// Following were commented out in python code. Let's keep them in, just in case.
	// Probably will return when we look at RNAPZ6 (adenosylcobalamin riboswitch),
	//  where Mg(2+) served as a virtual atom placeholder for adenosylcobalamin
	/////////////////////////////////////////////////////////////////////////////////////
	// need to handle Mg(2+)
	//mg_pos = []
	//for i in range( len( sequence ) ):
	//    if sequence[i]=='z': mg_pos.append( i+1 )
	//working_mg_pos = working_res_map( mg_pos, working_res )
	//if len( working_mg_pos ) > 0:
	//    for i in mg_pos:
	//        cutpoint_open_res_chain.append( i-1 )
	//        virtual_anchor.append( i )


	////////////////////
	// Step 16B
	////////////////////
	RNA_BasePairList rna_pairing_list;
	utility::vector1 < utility::vector1 <core::Size > > working_obligate_pair_sets;
	utility::vector1 < utility::vector1 <core::Size > > working_stem_pairing_sets;
	Size count( 0 );
	for ( Size n = 1; n <= working_stems.size(); n++ ) {
		vector1< std::pair< Size, Size > > working_stem = working_stems[ n ];
		vector1< Size > stem_pairing_set;
		for ( Size m = 1; m <= working_stem.size(); m++ ) {
			count++;
			BasePair base_pair( working_stem[m].first, working_stem[m].second, WATSON_CRICK, WATSON_CRICK, ANTIPARALLEL );
			runtime_assert( !rna_pairing_list.has_value( base_pair ) );
			rna_pairing_list.push_back( base_pair );
			stem_pairing_set.push_back( count );
		}
		working_stem_pairing_sets.push_back( stem_pairing_set );
	}

	for ( Size n = 1; n <= working_obligate_pairs.size(); n++ ) {
		BasePair const & base_pair = working_obligate_pairs[ n ];
		count++;
		rna_pairing_list.push_back( base_pair );
		working_obligate_pair_sets.push_back( make_vector1( count ) );
	}

	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;


	////////////////////
	// Step 17
	////////////////////
	///////////////////////////////////
	// package above into params
	///////////////////////////////////
	rna_params_ = RNA_DeNovoParametersOP( new RNA_DeNovoParameters );
	rna_params_->set_rna_pairing_list( rna_pairing_list );
	rna_params_->set_obligate_pairing_sets( working_obligate_pair_sets );
	rna_params_->set_stem_pairing_sets( working_stem_pairing_sets );
	rna_params_->set_chain_connections( working_chain_connections );
	rna_params_->set_cutpoints_open( working_cutpoint_open );
	rna_params_->set_cutpoints_closed( working_cutpoint_closed );
	rna_params_->set_virtual_anchor_attachment_points( working_virtual_anchor );


	////////////////////
	// Step 18
	////////////////////
	////////////////////////
	// Working native pose
	////////////////////////
	//Read in native if it exists.
	if ( option[ in::file::native ].user() ) {
		//Read in native if it exists.
		std::string native_pdb_file  = option[ in::file::native ]();
		native_pose_ = PoseOP( new Pose );
		core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
		pdbslice( *native_pose_, working_res );
	} else if ( option[ OptionKeys::rna::farna::working_native ].user() ) {
		std::string native_pdb_file  = option[ OptionKeys::rna::farna::working_native ];
		native_pose_ = PoseOP( new Pose );
		core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
	} else {
		runtime_assert( !option[ OptionKeys::rna::farna::refine_native ]() );
	}

	runtime_assert( option[ OptionKeys::rna::farna::minimize_rna ].user() ); // user should specify -minimize_rna true or -minimize_rna false
	// runtime_assert( option[ OptionKeys::score::include_neighbor_base_stacks ].user() ); // user should specify -include_neighbor_base_stacks true or -include_neighbor_base_stacks false.

	// some stuff to update in *options*
	options_->set_input_res( working_input_res );
	options_->set_extra_minimize_res( working_res_map( extra_minimize_res, working_res ) );
	options_->set_extra_minimize_chi_res( working_res_map( extra_minimize_chi_res, working_res ) );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// Following is copy of what was originally in rna_denovo.cc (main application)
// Should deprecate this in 2018, along with -use_legacy_setup option.
/////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoSetup::de_novo_setup_from_command_line_legacy()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::pose;
	using namespace protocols::farna::options;

	if ( options_->rna_params_file().size() > 0 )  {
		TR << TR.Red << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
		TR << TR.Red << " -params_file input will be deprecated soon. Use -secstruct_file, -obligate_pair, etc." << std::endl;
		TR << TR.Red << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
	}

	std::string const in_path = option[ in::path::path ]()[1];

	if ( option[ in::file::native ].user() ) {
		//Read in native if it exists.
		native_pose_ = pose::PoseOP( new pose::Pose );
		pose::Pose & native_pose = *native_pose_;
		std::string native_pdb_file  = option[ in::file::native ]();
		core::import_pose::pose_from_file( native_pose, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
	} else {
		runtime_assert( !option[ OptionKeys::rna::farna::refine_native ]() );
	}

	//Prepare starting structure from scratch --> read from fasta.
	pose_ = PoseOP( new Pose );
	if ( option[ OptionKeys::rna::farna::refine_native ]() ) {
		pose_ = native_pose_->clone();
	} else {
		std::string const sequence = core::sequence::read_fasta_file_return_str( option[ in::file::fasta ]()[1] );
		core::pose::make_pose_from_sequence( *pose_, sequence, *rsd_set_ );
	}

	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[ OptionKeys::constraints::cst_file ](1), ConstraintSetOP( new ConstraintSet ), *pose_ );
		pose_->constraint_set( cst_set );
	}

	if ( option[ OptionKeys::rna::data_file].user() ) {
		core::io::rna::RNA_DataReader rna_data_reader( in_path + option[ OptionKeys::rna::data_file ]  );
		rna_data_reader.fill_rna_data_info( *pose_ );
	}

	rna_params_ = RNA_DeNovoParametersOP( new RNA_DeNovoParameters( options_->rna_params_file() ) );

}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Utility functions
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void
RNA_DeNovoSetup::setup_refine_pose_list() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;

	if ( option[ OptionKeys::rna::farna::refine_native ]() ) options_->set_refine_pose( true );

	// Silent file input for fine refinement
	refine_pose_list_ = get_refine_pose_list(
		option[ OptionKeys::rna::farna::refine_silent_file ](),
		option[ OptionKeys::rna::farna::output_res_num ].resnum_and_chain(),
		rsd_set_ );

	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set( pose_->constraint_set()->clone() ); // assume constraints have been set up already
		for ( Size i = 1; i <= refine_pose_list_.size(); ++i ) refine_pose_list_[ i ]->constraint_set( cst_set );
	}

}
vector1<pose::PoseOP>
RNA_DeNovoSetup::get_refine_pose_list( std::string const & input_silent_file,
	std::pair< utility::vector1< int >, utility::vector1< char > > const & output_res_and_chain,
	core::chemical::ResidueTypeSetCOP rsd_set_ ) const
{
	vector1<pose::PoseOP> refine_pose_list;
	if ( input_silent_file.size() > 0 ) {
		core::import_pose::pose_stream::SilentFilePoseInputStream input( input_silent_file );
		input.set_order_by_energy( true );
		while ( input.has_another_pose() ) {
			pose::PoseOP new_pose( new pose::Pose );
			input.fill_pose( *new_pose, *rsd_set_ );
			protocols::farna::set_output_res_and_chain( *new_pose, output_res_and_chain );
			refine_pose_list.push_back( new_pose );
		}
	}
	return refine_pose_list;
}


vector1< Size >
RNA_DeNovoSetup::working_res_map( vector1< Size > const & vec,
	vector1< Size > const & working_res ) const
{
	if ( working_res.size() == 0 ) return vec;
	vector1< Size > working_vec;
	for ( Size i = 1; i <= vec.size(); i++ ) {
		Size const m = vec[ i ];
		if ( working_res.has_value( m ) ) {
			working_vec.push_back( working_res.index( m ) );
		}
	}
	return working_vec;
}

std::string
RNA_DeNovoSetup::working_res_map( std::string const & seq_input,
	vector1< Size > const & working_res ) const
{
	std::string seq( seq_input );
	core::sequence::strip_spacers( seq );
	if ( working_res.size() == 0 ) return seq;
	std::string working_seq;
	for ( Size m = 1; m <= seq.size(); m++ ) {
		if ( working_res.has_value( m ) ) {
			working_seq += seq[ m - 1 ];
		}
	}
	return working_seq;
}

// Following not handling spacers correctly...
secstruct::RNA_SecStruct
RNA_DeNovoSetup::working_res_map( secstruct::RNA_SecStruct const & rna_secstruct,
	vector1< Size > const & working_res ) const
{
	std::string working_secstruct = working_res_map( rna_secstruct.secstruct(), working_res );
	return secstruct::RNA_SecStruct( working_secstruct );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add by FCC: process PDBs and generate reasonable obligate pairs
void
RNA_DeNovoSetup::get_seq_and_resnum( std::string const & pdb,
	std::string & seq,
	vector1< int > & resnum,
	vector1< char > & chain ) const
{
	using namespace core::pose;
	using namespace core::import_pose;
	PoseOP pose_op = pose_from_file( pdb );
	Pose & pose = *pose_op;
	PDBInfoOP pdb_info = pose.pdb_info();
	seq = pose.sequence();
	resnum.clear();
	chain.clear();
	for ( Size n = 1; n <= pose.size(); n++ ) {
		resnum.push_back( pdb_info->number( n ) );
		chain.push_back( pdb_info->chain( n ) );
	}
}

std::string
RNA_DeNovoSetup::get_silent_seq( std::string const & silent_file ) const
{
	using namespace core::io::silent;
	SilentFileData silent_file_data;
	return silent_file_data.get_sequence( silent_file );
}

bool
RNA_DeNovoSetup::already_listed_in_obligate_pair( vector1< Size > const & new_pair,
	vector1< Size > const & obligate_pair ) const
{
	Size const new_pos1 = new_pair[ 1 ];
	Size const new_pos2 = new_pair[ 2 ];
	for ( Size m = 0; m < obligate_pair.size()/2; m++ ) {
		Size pos1 = obligate_pair[ 2*m+1 ];
		Size pos2 = obligate_pair[ 2*m+2 ];
		if ( pos1 == new_pos1 && pos2 == new_pos2 ) {
			return true;
		}
	}
	return false;
}

bool
RNA_DeNovoSetup::already_listed_in_obligate_pair( vector1< Size > const & new_pair,
	vector1< Size > const & obligate_pair,
	vector1< Size > const & obligate_pair_explicit ) const
{
	vector1< Size > all_pair;
	for ( Size n = 1; n <= obligate_pair.size(); n++ ) all_pair.push_back( obligate_pair[ n ] );
	for ( Size n = 1; n <= obligate_pair_explicit.size(); n++ ) all_pair.push_back( obligate_pair_explicit[ n ] );
	return already_listed_in_obligate_pair( new_pair, all_pair );
}

///////////////////////////////////////////////////////////////////
void
RNA_DeNovoSetup::update_working_obligate_pairs_with_stems(
	vector1< pose::rna::BasePair > & working_obligate_pairs,
	vector1< vector1< std::pair< Size, Size > > > const & working_stems,
	vector1< Size > const & working_input_res ) const
{
	for ( Size n = 1; n <= working_stems.size(); n++ ) {
		vector1< std::pair< Size, Size > > const & working_stem = working_stems[ n ];

		bool stem_in_input_res( false );
		for ( Size i = 1; i <= working_stem.size(); i++ ) {
			std::pair< Size, Size > const & pair = working_stem[ i ];
			if ( working_input_res.has_value( pair.first ) &&
					working_input_res.has_value( pair.second ) ) {
				stem_in_input_res = true; break;
			}
		}
		if ( stem_in_input_res ) continue;

		core::pose::rna::BasePair base_pair( working_stem[1].first, working_stem[1].second, WATSON_CRICK, WATSON_CRICK, ANTIPARALLEL );
		working_obligate_pairs.push_back( base_pair );
	}
}

///////////////////////////////////////////////////////////////////
vector1< std::pair< Size, Size > >
RNA_DeNovoSetup::flatten( vector1< vector1< std::pair< Size, Size > > > const & vec) const
{
	vector1< std::pair< Size, Size > > new_vec;
	for ( Size n = 1; n <= vec.size(); n++ ) {
		for ( Size m = 1; m <= vec[ n ].size(); m++ ) {
			new_vec.push_back( vec[ n ][ m ] );
		}
	}
	return new_vec;
}


} //setup
} //farna
} //protocols
