// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/setup/RNA_DeNovoSetup.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/denovo/setup/RNA_DeNovoSetup.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoParameters.hh>
#include <protocols/rna/denovo/options/RNA_DeNovoProtocolOptions.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/rna/movers/RNA_HelixAssembler.hh>
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
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <core/pose/rna/RNA_SecStruct.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <utility/options/OptionCollection.hh>

#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.denovo.setup.RNA_DeNovoSetup" );
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
namespace rna {
namespace denovo {
namespace setup {


using namespace core::pose::rna;
void
dump_stems( RNA_SecStruct const & working_secstruct, std::string const & working_sequence, core::pose::full_model_info::FullModelParameters const & full_model_parameters ) {

	movers::RNA_HelixAssembler rha;

	vector1< vector1< std::pair< Size, Size > > > const & working_stems = working_secstruct.stems();

	Size ii = 1;
	for ( auto const & stem : working_stems ) {
		std::stringstream filename_stream, helix_seq_stream;
		filename_stream << "helix_" << ii << ".pdb";
		utility::vector1< int > resnums;
		utility::vector1< char > chains;
		utility::vector1< std::string > segids;
		for ( Size ii = 1; ii <= stem.size(); ++ii ) {
			helix_seq_stream << working_sequence[stem[ii].first-1];
			auto res_chain_segid = full_model_parameters.full_to_conventional_resnum_and_chain_and_segid( stem[ii].first );
			resnums.push_back( std::get< 0 >( res_chain_segid ) );
			chains.push_back( std::get< 1 >( res_chain_segid ) );
			segids.push_back( std::get< 2 >( res_chain_segid ) );
		}
		for ( Size ii = stem.size(); ii >= 1; --ii ) {
			helix_seq_stream << working_sequence[stem[ii].second-1];
			auto res_chain_segid = full_model_parameters.full_to_conventional_resnum_and_chain_and_segid( stem[ii].second );
			resnums.push_back( std::get< 0 >( res_chain_segid ) );
			chains.push_back( std::get< 1 >( res_chain_segid ) );
			segids.push_back( std::get< 2 >( res_chain_segid ) );
		}

		pose::Pose helix_pose;
		rha.apply( helix_pose, helix_seq_stream.str() );

		pose::PDBInfoOP pdb_info( new pose::PDBInfo( helix_pose ) );
		pdb_info->set_chains( chains );
		pdb_info->set_numbering( resnums );
		pdb_info->set_segment_ids( segids );
		helix_pose.pdb_info( pdb_info );

		helix_pose.dump_pdb( filename_stream.str() );
		++ii;
	}
}


//Constructor
RNA_DeNovoSetup::RNA_DeNovoSetup():
	options_( options::RNA_DeNovoProtocolOptionsOP( new options::RNA_DeNovoProtocolOptions ) )
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
	using namespace protocols::rna::denovo::options;

	options_->initialize_from_command_line();

	rsd_set_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	if ( options_->use_legacy_setup() || options_->rna_params_file().size() > 0 ) {
		de_novo_setup_from_command_line_legacy();
	} else {
		de_novo_setup_from_command_line();
	}

	// if output_res_num supplied, this will change PDBInfo numbering & chain.
	set_output_res_and_chain( *pose_, option[ OptionKeys::rna::denovo::output_res_num ].resnum_and_chain() );

	// refine_pose is a seldom-used functionality at the moment -- not well tested.
	setup_refine_pose_list( option );

}

void
RNA_DeNovoSetup::initialize_from_options( utility::options::OptionCollection const & opts )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace protocols::rna::denovo::options;

	options_->initialize_from_options( opts );

	rsd_set_ = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	de_novo_setup_from_options( opts );

	// if output_res_num supplied, this will change PDBInfo numbering & chain.
	set_output_res_and_chain( *pose_, opts[ OptionKeys::rna::denovo::output_res_num ].resnum_and_chain() );

	// refine_pose is a seldom-used functionality at the moment -- not well tested.
	setup_refine_pose_list( opts );
}

void RNA_DeNovoSetup::initialize_inputs_from_options( utility::options::OptionCollection const & opts ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	fasta_files_ = opts[ in::file::fasta ]();
	input_pdbs_ = opts[ in::file::s ]();
	input_silent_files_ = opts[ in::file::silent ]();
	sequence_strings_ = opts[ OptionKeys::rna::denovo::sequence ]();
	if ( opts[ OptionKeys::rna::denovo::minimize_rna ].user() ) {
		minimize_rna_ = opts[ OptionKeys::rna::denovo::minimize_rna ]();
		minimize_rna_has_been_specified_ = true;
	}

}

void RNA_DeNovoSetup::set_input_pdbs( utility::vector1< std::string > const & input_pdbs ) {
	input_pdbs_ = input_pdbs;
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_chunk_pdb_files( input_pdbs_ );
}

void RNA_DeNovoSetup::set_input_silent_files( utility::vector1< std::string > const & input_silent_files ) {
	input_silent_files_ = input_silent_files;
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_chunk_silent_files( input_silent_files_ );
}

void RNA_DeNovoSetup::set_align_pdb( std::string const & align_pdb ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_align_pdb( align_pdb );
}

void RNA_DeNovoSetup::set_nstruct( Size const nstruct ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_nstruct( nstruct );
}

void RNA_DeNovoSetup::set_silent_file( std::string const & silent_file ) {
	runtime_assert( options_ );
	// This must propagate to options.
	options_->set_silent_file( silent_file );
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
	de_novo_setup_from_options( basic::options::option );
}

void
RNA_DeNovoSetup::de_novo_setup_from_options( utility::options::OptionCollection const & opts )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::rna::denovo::options;
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
	int const offset = opts[ OptionKeys::rna::denovo::offset ]();

	// Sequence setup:
	// FullModelParameters is a nice object that holds sequence, non-standard residues ("Z[Mg]"),
	// and what chains and residue numbers to use.
	FullModelParametersOP full_model_parameters;
	vector1< Size > cutpoint_open_in_full_model;
	if ( !fasta_files_.empty() ) {
		// use fasta readin developed for stepwise application -- also reads in
		// numbers & chains based on fasta header lines.
		if ( sequence_strings_.size() > 0 ) utility_exit_with_message( "Cannot specify both -sequence and -fasta" );
		if ( fasta_files_.size() != 1 ) utility_exit_with_message( "Please specify exactly one fasta file." );
		full_model_parameters = core::import_pose::get_sequence_information( fasta_files_[ 1 ], cutpoint_open_in_full_model );
		if ( offset != 0 ) {
			vector1< int > new_numbering = full_model_parameters->conventional_numbering();
			for ( Size n = 1; n <= new_numbering.size(); n++ ) { new_numbering[ n ] += offset; }
			full_model_parameters->set_conventional_numbering( new_numbering );
		}

	} else {
		// basic read-in of sequence from command line
		if ( sequence_strings_.size() == 0 ) utility_exit_with_message( "Must specify -sequence or -fasta" );
		std::string sequence( sequence_strings_[1] );
		for ( Size n = 2; n <= sequence_strings_.size(); n++ ) sequence += std::string( " " + sequence_strings_[ n ] );
		cutpoint_open_in_full_model = core::sequence::strip_spacers( sequence );
		std::map< Size, std::string > non_standard_residue_map = core::sequence::parse_out_non_standard_residues( sequence );
		vector1< int > res_numbers_in_pose;
		for ( Size n = 1; n <= sequence.size(); n++ ) res_numbers_in_pose.push_back( n + offset );
		core::import_pose::get_extra_cutpoints_from_names( sequence.size(), cutpoint_open_in_full_model, non_standard_residue_map );
		full_model_parameters = FullModelParametersOP( new FullModelParameters( sequence, cutpoint_open_in_full_model, res_numbers_in_pose ) );
		full_model_parameters->set_non_standard_residue_map( non_standard_residue_map );
		Size chain_num( 1 ), res( 0 );
		utility::vector1< char > chains; utility::vector1< Size > resnum;
		for ( Size n = 1; n <= sequence.size(); n++ ) {
			chains.push_back( chr_chains[ (chain_num - 1)  % chr_chains.size() ] );
			resnum.push_back( ++res );
			if ( cutpoint_open_in_full_model.has_value( n ) ) {
				chain_num++; res = 0;
			}
		}
		full_model_parameters->set_conventional_chains( chains );
		full_model_parameters->set_conventional_numbering( resnum );
	}
	std::string const sequence = full_model_parameters->full_sequence();

	////////////////////
	// Step 2
	////////////////////
	// Other useful residues.
	if ( opts[ full_model::cutpoint_open ].user() ) {
		cutpoint_open_in_full_model  =
			full_model_parameters->conventional_to_full( opts[ full_model::cutpoint_open ].resnum_and_chain() );
	}
	vector1< Size > working_res        =
		full_model_parameters->conventional_to_full( opts[ full_model::working_res ].resnum_and_chain() ); //all working stuff
	std::sort( working_res.begin(), working_res.end() ); // some the following depends on correct order.
	vector1< Size > const cutpoint_closed          =
		full_model_parameters->conventional_to_full( opts[ full_model::cutpoint_closed ].resnum_and_chain() );
	vector1< Size > const cutpoint_cyclize          =
		full_model_parameters->conventional_to_full( option[ full_model::cyclize ].resnum_and_chain() );
	// Ends up as pairs. Starts as a vector
	vector1< Size > const fiveprime_cap =
		full_model_parameters->conventional_to_full( option[ full_model::fiveprime_cap ].resnum_and_chain() );
	vector1< Size > block_stack_above_res  =
		full_model_parameters->conventional_to_full( opts[ full_model::rna::block_stack_above_res ].resnum_and_chain() );
	vector1< Size > block_stack_below_res  =
		full_model_parameters->conventional_to_full( opts[ full_model::rna::block_stack_below_res ].resnum_and_chain() );
	vector1< Size > extra_minimize_res =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::minimize::extra_minimize_res ].resnum_and_chain() );
	vector1< Size > extra_minimize_chi_res =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::minimize::extra_minimize_chi_res ].resnum_and_chain() );
	vector1< Size > input_res_user_defined =
		full_model_parameters->conventional_to_full( opts[ in::file::input_res ].resnum_and_chain() );
	vector1< Size > input_silent_res_user_defined =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::input_silent_res ].resnum_and_chain() );
	vector1< Size > virtual_anchor =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::virtual_anchor ].resnum_and_chain() );
	vector1< Size > obligate_pair =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::obligate_pair ].resnum_and_chain() );
	vector1< Size > remove_pair =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::remove_pair ].resnum_and_chain() );
	vector1< Size > remove_obligate_pair =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::remove_obligate_pair ].resnum_and_chain() );
	vector1< Size > output_jump_res =
		full_model_parameters->conventional_to_full( opts[ OptionKeys::rna::denovo::out::output_jump_res ].resnum_and_chain() );

	////////////////////
	// Step 3
	////////////////////
	// secondary structure setup.
	RNA_SecStruct secstruct( opts[ OptionKeys::rna::denovo::secstruct ](), opts[ OptionKeys::rna::denovo::secstruct_file ](), sequence );
	// "general" secondary structure includes non-canonical pairs that should be connected by jumps during run; used with -bps_moves.
	RNA_SecStruct secstruct_general( opts[ OptionKeys::rna::denovo::secstruct_general ](), opts[ OptionKeys::rna::denovo::secstruct_general_file ](), sequence );

	secstruct.check_compatible_with_sequence( sequence, true  /*check_complementarity*/ );
	secstruct_general.check_compatible_with_sequence( sequence, false /*check_complementarity*/ );

	if ( !secstruct_general.blank() && !options_->bps_moves() ) utility_exit_with_message("cannot supply secstruct_general without bps_moves");

	////////////////////
	// Step 4
	////////////////////
	/////////////
	std::string const working_native_pdb = opts[ OptionKeys::rna::denovo::working_native ]();
	vector1< std::string > obligate_pair_explicit = opts[ OptionKeys::rna::denovo::obligate_pair_explicit ]();
	vector1< std::string > chain_connections = opts[ OptionKeys::rna::denovo::chain_connection ]();

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
	for ( std::string const & pdb : input_pdbs_ ) {
		std::string pdb_seq;
		vector1< int > resnum;
		vector1< char >  chain;
		vector1< std::string >  segid;
		get_seq_and_resnum( pdb, pdb_seq, resnum, chain, segid );
		vector1< Size > resnum_in_full_model;

		if ( input_res_user_defined_count + resnum.size() <= input_res_user_defined.size() ) {
			// input res could have come from user after flag -input_res
			for ( Size q = 1; q <= resnum.size(); q++ ) {
				input_res_user_defined_count++;
				resnum_in_full_model.push_back( input_res_user_defined[ input_res_user_defined_count ] );
			}
		} else {
			// figure out residue numbers from PDB resnum & chain.
			resnum_in_full_model = full_model_parameters->conventional_to_full( std::make_tuple( resnum, chain, segid ) );
		}

		std::string actual_seq = "";
		for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
			Size const i = resnum_in_full_model[ q ];
			if ( input_res.has_value( i ) )  TR << TR.Red << "WARNING! Input residue " << resnum[q] << " " << chain[q] << " exists in two pdb files!!" << std::endl;
			actual_seq += sequence[ i - 1 ];
			input_res.push_back( i );
		}
		if ( pdb_seq != actual_seq ) {
			TR << TR.Red << "   pdb_seq: " << pdb_seq << std::endl;
			TR << TR.Red << "target_seq: " << actual_seq << std::endl;
			for ( Size q = 1; q <= resnum_in_full_model.size(); q++ ) {
				if ( q > actual_seq.size() ) {
					TR << "mismatch in length beyond " << q << std::endl;
					break;
				}
				if ( pdb_seq[q-1] != actual_seq[q-1] ) {
					Size n( resnum_in_full_model[ q ] );
					TR << "mismatch in sequence: pdb  " << pdb_seq[q-1] << " vs target " << actual_seq[q-1] << " at " << full_model_parameters->full_to_conventional( n ) << std::endl;
				}
			}
			utility_exit_with_message("The sequence in "+pdb+" does not match target sequence!!");
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

	for ( std::string const & silent : input_silent_files_ ) {
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
			vector1<Size> input_silent_res_from_file = full_model_parameters->conventional_to_full( get_silent_resnum( silent ) );
			for ( Size const res : input_silent_res_from_file ) {
				resnum_in_full_model.push_back( res );
			}
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
		std::vector< std::string > segids;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1< std::string >( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1< std::string >( segids ) ) );
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
			if ( j > 0 &&
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

			Size segment1_end   = chunks[ i     ][ chunks[ i ].size() ];
			Size segment2_start = chunks[ i_next ][ 1 ];

			// avoid placing jump positions at extra_minimize_res
			while ( extra_minimize_res.has_value( segment1_end )   && segment1_end > chunks[ i ][ 1 ] ) segment1_end--;
			while ( extra_minimize_res.has_value( segment2_start ) && segment2_start < chunks[ i_next ][ chunks[ i_next ].size() ] ) segment2_start++;

			Size const new_pos1 = std::min( segment1_end, segment2_start );
			Size const new_pos2 = std::max( segment1_end, segment2_start );
			obligate_pair.push_back( new_pos1 );
			obligate_pair.push_back( new_pos2 );
			//   TR << TR.Cyan << "Creating new obligate pair: " << obligate_pair << " for chunk with residues " << resnum << std::endl;
			n_jumps++;
		}
	}
	vector1< std::pair< Size, Size > > canonical_pairs = secstruct.base_pairs();
	vector1< std::pair< Size, Size > > general_pairs   = secstruct_general.base_pairs();
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
	RNA_SecStruct working_secstruct = working_res_map( secstruct, working_res );
	RNA_SecStruct working_secstruct_general = working_res_map( secstruct_general, working_res );

	TR << "Sequence:            " << working_sequence << std::endl;
	TR << "Secstruct:           " << working_secstruct.secstruct() << std::endl;
	if ( !secstruct_general.blank() ) TR << "Secstruct [general]: " << working_secstruct_general.secstruct() << std::endl;

	// Step 8B [good ol' secstruct_legacy]
	std::string secstruct_legacy = opts[ OptionKeys::rna::denovo::secstruct_legacy ]();
	std::string working_secstruct_legacy = working_res_map( secstruct_legacy, working_res );
	if ( secstruct_legacy.size() > 0  ) TR << "Secstruct [legacy]: " << working_secstruct_legacy << std::endl;

	if ( options_->dump_stems() ) dump_stems( working_secstruct, working_sequence, *full_model_parameters );

	////////////////////
	// Step 9
	////////////////////
	///////////////////////////////////////////////////////////////////
	// initialize variables needed for RNA_DeNovoParams (rna_params_)
	///////////////////////////////////////////////////////////////////
	vector1< Size > working_cutpoint_open   = working_res_map( cutpoint_open_in_full_model, working_res, true /*leave out last residue*/ );
	vector1< Size > working_cutpoint_closed = working_res_map( cutpoint_closed, working_res );
	vector1< Size > working_cutpoint_cyclize = working_res_map( cutpoint_cyclize, working_res );
	vector1< Size > working_fiveprime_cap = working_res_map( fiveprime_cap, working_res );
	/*vector1< std::pair< Size, Size > > working_fiveprime_cap_in_pairs;
	for ( Size ii = 1; ii <= working_fiveprime_cap.size() - 1; ii += 2 ) {
	working_fiveprime_cap_in_pairs.emplace_back( std::make_pair( working_fiveprime_cap[ii], working_fiveprime_cap[ii+1] ) );
	}*/
	vector1< Size > working_block_stack_above_res = working_res_map( block_stack_above_res, working_res );
	vector1< Size > working_block_stack_below_res = working_res_map( block_stack_below_res, working_res );
	vector1< Size > working_virtual_anchor  = working_res_map( virtual_anchor, working_res );
	vector1< pose::rna::BasePair > working_obligate_pairs;
	vector1< vector1< std::pair< Size, Size > > > working_stems;

	for ( Size n = 1; n < working_res.size(); n++ ) {
		if ( working_res[ n+1 ] > working_res[n] + 1  && !working_cutpoint_open.has_value( n ) )  working_cutpoint_open.push_back( n );
	}

	/////////////////////////
	// working stems
	/////////////////////////
	vector1< Size > working_cutpoint( working_cutpoint_open ); working_cutpoint.append( working_cutpoint_closed );
	working_stems = working_secstruct.stems();
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
	std::string const in_path = opts[ in::path::path ]()[1];
	Pose full_pose;
	pose_ = PoseOP( new Pose );
	std::string const full_annotated_sequence = full_model_parameters->full_annotated_sequence();
	make_pose_from_sequence( full_pose, full_annotated_sequence, *rsd_set_ );
	set_output_res_and_chain( full_pose, std::make_tuple( full_model_parameters->conventional_numbering(),
		full_model_parameters->conventional_chains(), full_model_parameters->conventional_segids()  ) );
	// Check whether the sequence contains protein residues
	bool is_rna_and_protein = false;
	for ( core::Size r=1; r<=full_pose.total_residue(); ++r ) {
		if ( full_pose.residue( r ).is_protein() ) {
			is_rna_and_protein = true;
			break;
		}
	}

	//
	if ( is_rna_and_protein ) {
		if ( !opts[ OptionKeys::rna::denovo::lores_scorefxn ].user() ) {
			// set default low-res RNA/protein score function
			options_->set_lores_scorefxn( "rna/denovo/rna_lores_with_rnp" );
		}
	}

	////////////////////////
	// Working native pose
	////////////////////////
	//Read in native if it exists.
	if ( opts[ in::file::native ].user() ) {
		// AMW TODO: later, let this pass align_pose on to the RNA_DeNovoProtocol
		// and the RNA_FragmentMonteCarlo. At the moment, this isn't necessary at
		// all, though.
		PoseOP align_pose;
		stepwise::setup::initialize_native_and_align_pose( native_pose_, align_pose, rsd_set_, pose_ );

		// if refine native, set the full_pose equal to the native_pose
		if ( opts[ OptionKeys::rna::denovo::refine_native ]() ) {
			full_pose = *(native_pose_->clone());
		}
		// AMW: this function will ensure that it just copies the native pose if in fact working_res is 'everything'
		pdbslice( *native_pose_, working_res );
	} else if ( opts[ OptionKeys::rna::denovo::working_native ].user() ) {
		std::string native_pdb_file  = opts[ OptionKeys::rna::denovo::working_native ];
		native_pose_ = PoseOP( new Pose );
		core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
	} else {
		runtime_assert( !opts[ OptionKeys::rna::denovo::refine_native ]() );
	}

	// if the refine_native, then the full pose is equal to the native pose
	if ( ! opts[ OptionKeys::rna::denovo::working_native ].user() ) { // usually not defined by user
		pdbslice( *pose_, full_pose, working_res );
	} else {
		// there might still be issues with how the csts are set up here...?
		pose_ = native_pose_->clone();
	}



	////////////////////
	// Step 11
	////////////////////
	////////////////////////
	// working constraints.
	////////////////////////
	if ( opts[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set( new ConstraintSet );
		//opts[ OptionKeys::constraints::force_pdb_info_mapping ].def( true ); // using option as global variable due to difficulty in dealing with static functions.
		cst_set = ConstraintIO::get_instance()->read_constraints( opts[ OptionKeys::constraints::cst_file ](1), ConstraintSetOP( new ConstraintSet ), full_pose, true /*force_pdb_info_mapping*/ );
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
	if ( opts[ OptionKeys::rna::data_file].user() ) {
		core::io::rna::RNA_DataReader rna_data_reader( in_path + opts[ OptionKeys::rna::data_file ]  );
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
		std::vector< std::string > segids;
		vector1< Size > resnum_full;

		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 1 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1<std::string>( segids ) ) );
		runtime_assert( resnum_full.size() == 1 );
		Size const pos1 = resnum_full[ 1 ];
		if ( !working_res.has_value( pos1 ) ) continue;

		resnum.clear(); chains.clear();
		utility::get_resnum_and_chain_from_one_tag( obligate_pair_explicit[ 5*m + 2 ], resnum, chains, segids );
		resnum_full = full_model_parameters->conventional_to_full( std::make_tuple( vector1<int>( resnum ),
			vector1<char>( chains ), vector1<std::string>( segids ) ) );
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
				std::vector< std::string > segids;
				bool ok = utility::get_resnum_and_chain_from_one_tag( chain_connections[ k ], resnum, chains, segids );
				runtime_assert( ok );
				vector1< Size > resnum_full = full_model_parameters->conventional_to_full(
					std::make_tuple( vector1<Size>( resnum ),
					vector1< char >( chains ),
					vector1< std::string >( segids ) ) );
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
	rna_params_->set_cutpoints_cyclize( working_cutpoint_cyclize );
	rna_params_->set_fiveprime_cap( working_fiveprime_cap );
	rna_params_->set_block_stack_above_res( working_block_stack_above_res );
	rna_params_->set_block_stack_below_res( working_block_stack_below_res );
	rna_params_->set_virtual_anchor_attachment_points( working_virtual_anchor );
	rna_params_->set_rna_and_protein( is_rna_and_protein );
	rna_params_->set_rna_secstruct_legacy( working_secstruct_legacy );


	////////////////////
	// Step 18
	////////////////////
	// Moved this up to restore refine_native functionality
	//////////////////////////
	//// Working native pose
	//////////////////////////
	////Read in native if it exists.
	//if ( opts[ in::file::native ].user() ) {
	// //Read in native if it exists.
	// std::string native_pdb_file  = opts[ in::file::native ]();
	// native_pose_ = PoseOP( new Pose );
	// core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
	// pdbslice( *native_pose_, working_res );
	// // set the pose equal to the native pose if user wants to refine_native
	// // does having this way down here mess anything up?
	// if ( opts[ OptionKeys::rna::denovo::refine_native ]() ) {
	//  pose_ = native_pose_->clone();
	// }
	//} else if ( opts[ OptionKeys::rna::denovo::working_native ].user() ) {
	// std::string native_pdb_file  = opts[ OptionKeys::rna::denovo::working_native ];
	// native_pose_ = PoseOP( new Pose );
	// core::import_pose::pose_from_file( *native_pose_, *rsd_set_, in_path + native_pdb_file , core::import_pose::PDB_file);
	//} else {
	// runtime_assert( !opts[ OptionKeys::rna::denovo::refine_native ]() );
	//}

	if ( !opts[ OptionKeys::rna::denovo::minimize_rna ].user() && !minimize_rna_has_been_specified_ ) utility_exit_with_message( "Please specify either '-minimize_rna true' or '-minimize_rna false'." );

	if ( minimize_rna_has_been_specified_ ) {
		options_->set_minimize_structure( minimize_rna_ );
	}

	// runtime_assert( opts[ OptionKeys::score::include_neighbor_base_stacks ].user() ); // user should specify -include_neighbor_base_stacks true or -include_neighbor_base_stacks false.

	// some stuff to update in *options*
	options_->set_input_res( working_input_res );
	options_->set_extra_minimize_res( working_res_map( extra_minimize_res, working_res ) );
	options_->set_extra_minimize_chi_res( working_res_map( extra_minimize_chi_res, working_res ) );
	options_->set_output_jump_res( working_res_map( output_jump_res, working_res ) );

	////////////////////
	// Step 19
	////////////////////
	using namespace core::pose::full_model_info;
	// could also set up other stuff inside full_model_parameters -- see protocols/stepise/FullModelInfoSetupFromCommandLine.cc
	// a better route, however, would be to *deprecate* this setup code, and instead use that stepwise code + build_full_model to
	// handle FARFAR setup. -- rhiju & amwatkins, dec. 2016.
	full_model_parameters->set_parameter_as_res_list( CUTPOINT_OPEN, cutpoint_open_in_full_model );
	full_model_parameters->set_parameter_as_res_list( RNA_SYN_CHI,
		full_model_parameters->conventional_to_full( opts[ full_model::rna::force_syn_chi_res_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_ANTI_CHI,
		full_model_parameters->conventional_to_full( opts[ full_model::rna::force_anti_chi_res_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_ABOVE, block_stack_above_res );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_BELOW, block_stack_below_res );
	full_model_parameters->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res );

	vector1< Size > dummy_domain_map( sequence.size(), 0 );
	full_model_parameters->set_parameter( INPUT_DOMAIN, dummy_domain_map /* domain_map */ );
	full_model_parameters->set_parameter( FIXED_DOMAIN, dummy_domain_map /*stepwise::setup::figure_out_fixed_domain_map( domain_map, extra_minimize_res ) */  );

	// Set up FullModelInfo (so we can use info stored there like SYN_CHI_RES)
	FullModelInfoOP full_model_info( new FullModelInfo( full_model_parameters ) );
	full_model_info->set_res_list( working_res );
	set_full_model_info( *pose_, full_model_info );
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
	using namespace protocols::rna::denovo::options;

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
		runtime_assert( !option[ OptionKeys::rna::denovo::refine_native ]() );
	}

	//Prepare starting structure from scratch --> read from fasta.
	pose_ = PoseOP( new Pose );
	if ( option[ OptionKeys::rna::denovo::refine_native ]() ) {
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
RNA_DeNovoSetup::setup_refine_pose_list( utility::options::OptionCollection const & opts ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;

	if ( opts[ OptionKeys::rna::denovo::refine_native ]() ) options_->set_refine_pose( true );

	// Silent file input for fine refinement
	refine_pose_list_ = get_refine_pose_list(
		opts[ OptionKeys::rna::denovo::refine_silent_file ](),
		opts[ OptionKeys::rna::denovo::output_res_num ].resnum_and_chain(),
		rsd_set_ );

	if ( opts[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set( pose_->constraint_set()->clone() ); // assume constraints have been set up already
		for ( Size i = 1; i <= refine_pose_list_.size(); ++i ) refine_pose_list_[ i ]->constraint_set( cst_set );
	}
}

vector1<pose::PoseOP>
RNA_DeNovoSetup::get_refine_pose_list( std::string const & input_silent_file,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & output_res_and_chain_and_segid,
	core::chemical::ResidueTypeSetCOP rsd_set_ ) const
{
	vector1<pose::PoseOP> refine_pose_list;
	if ( !input_silent_file.empty() ) {
		core::import_pose::pose_stream::SilentFilePoseInputStream input( input_silent_file );
		input.set_order_by_energy( true );
		while ( input.has_another_pose() ) {
			pose::PoseOP new_pose( new pose::Pose );
			input.fill_pose( *new_pose, *rsd_set_ );
			protocols::rna::denovo::set_output_res_and_chain( *new_pose, output_res_and_chain_and_segid );
			refine_pose_list.push_back( new_pose );
		}
	}
	return refine_pose_list;
}


vector1< Size >
RNA_DeNovoSetup::working_res_map( vector1< Size > const & vec,
	vector1< Size > const & working_res,
	bool const leave_out_last_working_residue /* = false */ ) const
{
	if ( working_res.size() == 0 ) return vec;
	vector1< Size > working_vec;
	for ( Size const m : vec ) {
		if ( leave_out_last_working_residue && m == working_res[ working_res.size() ]  ) continue;
		if ( working_res.has_value( m ) ) {
			working_vec.push_back( working_res.index( m ) );
		}
	}
	return working_vec;
}

std::string
RNA_DeNovoSetup::working_res_map( std::string const & seq_input,
	vector1< Size > const & working_res,
	bool const annotations_in_brackets /* = true */ ) const
{
	if ( seq_input.size() == 0 ) return "";
	std::string seq( seq_input );
	core::sequence::strip_spacers( seq, annotations_in_brackets );
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
core::pose::rna::RNA_SecStruct
RNA_DeNovoSetup::working_res_map( core::pose::rna::RNA_SecStruct const & rna_secstruct,
	vector1< Size > const & working_res ) const
{
	std::string working_secstruct = working_res_map( rna_secstruct.secstruct(), working_res, false /*annotations_in_brackets*/ );
	return core::pose::rna::RNA_SecStruct( working_secstruct );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add by FCC: process PDBs and generate reasonable obligate pairs
void
RNA_DeNovoSetup::get_seq_and_resnum( std::string const & pdb,
	std::string & seq,
	vector1< int > & resnum,
	vector1< char > & chain,
	vector1< std::string > & segid ) const
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
		segid.push_back( pdb_info->segmentID( n ) );
	}
}

std::string
RNA_DeNovoSetup::get_silent_seq( std::string const & silent_file ) const
{
	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );
	return silent_file_data.get_sequence( silent_file );
}

std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > >
RNA_DeNovoSetup::get_silent_resnum( std::string const & silent_file ) const
{
	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );
	return silent_file_data.get_resnum( silent_file );
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
	vector1< Size > all_pair = obligate_pair;
	//for ( Size n = 1; n <= obligate_pair.size(); n++ ) all_pair.push_back( obligate_pair[ n ] );
	for ( Size const exp : obligate_pair_explicit ) all_pair.push_back( exp );
	return already_listed_in_obligate_pair( new_pair, all_pair );
}

///////////////////////////////////////////////////////////////////
void
RNA_DeNovoSetup::update_working_obligate_pairs_with_stems(
	vector1< pose::rna::BasePair > & working_obligate_pairs,
	vector1< vector1< std::pair< Size, Size > > > const & working_stems,
	vector1< Size > const & working_input_res ) const
{
	for ( auto const & working_stem : working_stems ) {
		bool stem_in_input_res( false );

		for ( auto const & pair : working_stem ) {
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

} //setup
} //denovo
} //rna
} //protocols
