// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/util.cc
/// @brief small bundle of utilities for dealing with sequences
/// @author James Thompson
/// @author Sergey Lyskov

// Unit header
#include <core/sequence/util.hh>

// C/C++ headers
#include <algorithm>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/scoring/rms_util.hh>

// Package headers
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/DerivedSequenceMapping.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
namespace core {
namespace sequence {

using utility::vector1;
using std::string;

static basic::Tracer tr( "core.sequence" );

// any of these are OK as strand separators
vector1< char > spacers = utility::tools::make_vector1( '+','*',' ',',' );

void read_all_alignments(std::string const & format,
	utility::vector1<std::string> const & files,
	utility::vector1<SequenceAlignment> & alignments) {
	using std::string;
	using utility::vector1;

	for ( auto const & file : files ) {
		vector1<SequenceAlignment> current = read_aln(format, file);
		std::copy(current.begin(), current.end(), std::back_inserter(alignments));
	}
}

///////////////////////////////////////////////////////////////////////////////
///
/// @details if position i in seq1 is aligned with position j in seq2, mapping[ i ] == j
/// if position i in seq1 is unaligned, mapping[ i ] == 0
//
void
read_alignment_file(
	std::string const & filename,
	std::string & seq1,
	std::string & seq2,
	sequence::DerivedSequenceMapping & mapping // from numbering in sequence 1 to numbering in sequence 2
)
{
	std::string align1, align2;
	{ // parse the file
		std::ifstream data( filename.c_str() );
		std::string line;
		// 1st sequence
		getline( data,line );
		runtime_assert( line[0] == '>' );
		getline( data, align1 );
		// 2nd sequence
		getline( data, line );
		runtime_assert( line[0] == '>' );
		getline( data, align2 );
		data.close();
	}

	runtime_assert( align1.size() == align2.size() );

	seq1.clear();
	seq2.clear();
	mapping.clear();
	int pos1(0), pos2(0);
	for ( Size i=0; i< align1.size(); ++i ) {
		char const al1( align1[i] ), al2( align2[i] );
		bool const gap1( al1 == '.' || al1 == '-' );
		bool const gap2( al2 == '.' || al2 == '-' );
		if ( !gap2 ) {
			++pos2;
			seq2 += al2;
		}
		if ( !gap1 ) {
			++pos1;
			seq1 += al1;
			if ( !gap2 ) {
				mapping.push_back( pos2 );
			} else {
				mapping.push_back( 0 ); // unaligned
			}
		}
	}

	tr << "DerivedSequenceMapping align1: " << align1 << "\nalign2: " << align2 <<
		"\nseq1: " << seq1 << "\nseq2: " << seq2 << '\n';

	runtime_assert( mapping.size1() == seq1.size() );
	mapping.size2( seq2.size() ); // set sequence 2 size
} // read_alignment_file

vector1< string > read_fasta_file_str( std::string const & filename ) {
	vector1< SequenceOP > seqs = read_fasta_file( filename );
	vector1< string > seq_strings;
	for ( SequenceOP seq : seqs ) {
		seq_strings.push_back( seq->sequence() );
	}
	return seq_strings;
}

SequenceOP
get_sequence_object( std::string const & current_id,
	std::string const & current_sequence )
{
	std::string current_id_strip ( current_id );
	ObjexxFCL::strip_whitespace( current_id_strip );

	std::string current_sequence_strip( current_sequence );
	vector1< Size > spacer_positions = strip_spacers( current_sequence_strip );

	SequenceOP sequence( new Sequence( current_sequence_strip, current_id_strip ) );
	sequence->spacer_positions( spacer_positions );
	return sequence;
}

vector1< SequenceOP > read_fasta_file( std::string const & filename ) {
	vector1< SequenceOP > sequences;

	utility::io::izstream input( filename.c_str() );
	std::string line, current_sequence, current_id = "empty";

	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
		return sequences;
	}

	while ( getline( input, line ) ) {
		if ( line.substr(0,1) == ">" ) {
			if ( current_sequence != "" ) {
				sequences.push_back( get_sequence_object( current_id, current_sequence ) );
				current_sequence = "";
			}
			current_id = line.substr(1,line.size());
			continue;
		}
		current_sequence = current_sequence + ObjexxFCL::rstrip(line);
	}
	if ( current_sequence != "" ) {
		sequences.push_back( get_sequence_object( current_id, current_sequence ) );
	}

	return sequences;
} // read_fasta_file


std::string read_fasta_file_return_str( std::string const & filename ) {
	utility::vector1< std::string > sequences = read_fasta_file_str( filename );
	std::string full_sequence;
	for ( Size n = 1; n <= sequences.size(); n++ ) full_sequence += sequences[ n ];
	return full_sequence;
} // read_fasta_file


///  @brief read sequence from particular section of fasta file (comment starting with '> section'), terminate with failure if section not found
///         Note: section detection string is case insensitive
std::string read_fasta_file_section(std::string const & filename, std::string const & section_)
{
	string section(section_);
	std::transform(section.begin(), section.end(), section.begin(), ::tolower);

	vector1< SequenceOP > S = read_fasta_file(filename);

	for ( Size i=1; i<=S.size(); ++i ) {
		string id = S[i]->id();
		std::transform(id.begin(), id.end(), id.begin(), ::tolower);
		if ( utility::startswith(id, section) ) return S[i]->sequence();
	}

	utility_exit_with_message( "Error can't find section '" + section + "' in fasta file " + filename + "!" );
	return "";
}

///////////////////////////////////////////////////////////////////////////////////////
// looks for tab-delimited tags like 'chain:A' and 'res_num:5-20' in fasta IDs.
///////////////////////////////////////////////////////////////////////////////////////
void
get_conventional_chains_and_numbering( utility::vector1< SequenceCOP > const & fasta_sequences,
	utility::vector1< char > & conventional_chains,
	utility::vector1< int > & conventional_numbering,
	utility::vector1< std::string > & conventional_segids ) {
	using utility::string_split;
	bool found_info_in_previous_sequence( false );
	Size count( 0 );
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		utility::vector1< char > chains;
		utility::vector1< int > resnum;
		utility::vector1< std::string > segids;
		bool found_info( false );
		std::string tag;
		std::stringstream ss( fasta_sequences[n]->id() );
		while ( ss.good() ) {
			ss >> tag;
			bool string_is_ok( false );
			std::tuple< std::vector< int >, std::vector< char >, std::vector< std::string > > resnum_and_chain_and_segid = utility::get_resnum_and_chain_and_segid( tag, string_is_ok );
			if ( !string_is_ok ) continue;
			for ( Size n = 0; n < std::get< 0 >( resnum_and_chain_and_segid ).size(); n++ ) {
				if ( std::get< 1 >( resnum_and_chain_and_segid )[n] == ' ' ) continue; // there better be a chain. accept "A:1" but not just "1"
				resnum.push_back( std::get< 0 >( resnum_and_chain_and_segid )[n] );
				chains.push_back( std::get< 1 >( resnum_and_chain_and_segid )[n] );
				segids.push_back( std::get< 2 >( resnum_and_chain_and_segid )[n] );
			}
			found_info = true;
		}
		if ( n > 1 ) runtime_assert( found_info == found_info_in_previous_sequence );

		Size const clean_len = core::pose::rna::remove_bracketed( fasta_sequences[n]->sequence() ).size();
		if ( !found_info || resnum.size() != clean_len ) { /*happens with stray numbers*/
			resnum.clear();
			for ( Size q = 1; q <= clean_len; ++q ) {
				resnum.push_back( ++count );
				chains.push_back( ' ' ); // unknown chain
				segids.push_back( "    " ); // unknown chain
			}
		}
		std::string const sequence = fasta_sequences[n]->sequence();
		runtime_assert( clean_len == resnum.size() ); //sequence.size() == resnum.size() );
		for ( Size q = 1; q <= clean_len; q++ ) conventional_chains.push_back( chains[ q ] );
		for ( Size q = 1; q <= clean_len; q++ ) conventional_numbering.push_back( resnum[q] );
		for ( Size q = 1; q <= clean_len; q++ ) conventional_segids.push_back( segids[q] );

		found_info_in_previous_sequence  = found_info;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
/// @brief Return a string of concatenated SequenceCOP sequences
/// @details moved from stepwise/setup/FullModelInfoSetupFromCommandLine.cc
///////////////////////////////////////////////////////////////////////////////////////
std::string
get_concatenated_sequence( vector1< SequenceCOP > const & fasta_sequences ) {
	std::string sequence;
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		sequence += fasta_sequences[n]->sequence();
	}
	return sequence;
}

///////////////////////////////////////////////////////////////////////////////////////
/// @brief Read fasta file and concatenate sequences
///////////////////////////////////////////////////////////////////////////////////////
std::string
read_fasta_file_and_concatenate( std::string const & filename ) {
	vector1< SequenceOP > fasta_sequences = read_fasta_file( filename );
	std::string sequence = get_concatenated_sequence( fasta_sequences );
	return sequence;
}

core::sequence::DerivedSequenceMapping simple_mapping_from_file( std::string const & filename ) {
	// file I/O stuff
	utility::io::izstream input( filename.c_str() );
	std::string line;

	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
	}

	utility::vector1< std::pair< Size,Size > > aligned;
	Size max_resi = 0, max_resj = 0;

	std::string seq1, seq2;
	Size start_seq2( 0 );
	while ( getline( input, line ) ) {
		//can I read ungapped ?
		tr.Trace << "read line: " << line << std::endl;
		{
			std::istringstream line_stream( line );
			std::string tag;
			line_stream >> tag;
			if ( line_stream && tag.substr(0, std::string("ungapped").size() ) == "ungapped" ) {
				std::string type = tag.substr( std::string("ungapped_").size() );
				if ( type == "template:" ) {
					line_stream >> seq2;
					tr.Info << "read template sequence " << seq2 << std::endl;
				} else if ( type == "query:" ) {
					line_stream >> seq1;
					tr.Info << "read query sequence " << seq1 << std::endl;
				} else {
					utility_exit_with_message( "expected either ungapped_template or ungapped_query in file " + filename );
				}
				continue; // next line
			} // read sequence
		} //scope
		std::istringstream line_stream( line );
		Size resi, resj;
		line_stream >> resi >> resj;

		aligned.push_back( std::make_pair( resi, resj ) );
		max_resi = std::max( resi, max_resi );
		max_resj = std::max( resj, max_resj );
		if ( start_seq2 == 0 ) start_seq2 = resj;
	} // while( getline( input, line ) )

	// create SequenceMapping object
	DerivedSequenceMapping mapping( max_resi, max_resj );
	for ( vector1< std::pair< Size,Size > >::const_iterator it = aligned.begin(), end = aligned.end();
			it != end; ++it ) {
		mapping.insert_aligned_residue_safe( it->first, it->second );
	} // for ( aligned )

	mapping.seq1( seq1 );
	mapping.seq2( seq2 );
	mapping.start_seq2( start_seq2 );

	// if ( tr.Trace.visible() ) mapping.show( tr.Trace );

	return mapping;
} // mapping_from_file

utility::vector1< SequenceOP > seqs_from_cmd_lines() {
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< SequenceOP > seqs;

	if ( option[ in::file::fasta ].user() ) {
		vector1< FileName > fns( option[ in::file::fasta ]() );
		for ( FileName const & fn : fns ) {
			vector1< SequenceOP > temp_seqs( read_fasta_file( fn ) );
			for ( SequenceOP seq : temp_seqs ) {
				seqs.push_back( seq );
			}
		}
	}

	if ( option[ in::file::pssm ].user() ) {
		vector1< FileName > fns( option[ in::file::pssm ]() );
		for ( FileName const & fn : fns ) {
			SequenceProfileOP prof( new SequenceProfile );
			prof->read_from_file( fn );
			prof->convert_profile_to_probs(); // was previously implicit in read_from_file()
			seqs.push_back( prof );
		}
	}

	return seqs;
}

utility::vector1< SequenceAlignment > read_aln(
	std::string const & format,
	std::string const & filename
) {

	utility::vector1< SequenceAlignment > retval;
	if ( format == "general" ) {
		retval = read_general_aln_file( filename );
	} else if ( format == "grishin" ) {
		retval = read_grishin_aln_file( filename );
	} else {
		utility_exit_with_message(
			std::string( "No match for format " + format + "!" )
		);
	}
	tr.Debug << "read " << retval.size() << " alignments from file " << filename
		<< " with format " << format << "." << std::endl;
	return retval;
}


utility::vector1< SequenceAlignment > read_general_aln(
	std::istream & input
) {
	std::string line;
	utility::vector1< SequenceAlignment > alignments;
	SequenceAlignmentOP current( new SequenceAlignment );
	while ( getline( input, line ) ) {
		if ( line.substr(0,5) == "score" ) {
			std::istringstream line_input( line );
			std::string dummy;
			Real score( 0.0 );

			line_input >> dummy >> score;
			current->score( score );
			//while ( line_input.is_good() ) {
			//
			//}
		} else if ( line.substr(0,2) == "--" ) {
			if ( current->size() > 0 ) alignments.push_back( *current );
			current = SequenceAlignmentOP( new SequenceAlignment );
		} else if ( line.substr(0,1) == "#" ) {
			// do nothing. eventually add a comment here.
		} else {
			std::istringstream line_input( line );
			SequenceOP new_seq( new Sequence );
			new_seq->read_data( line_input );
			current->add_sequence( new_seq );
		}
	} // while
	if ( current->size() > 0 ) {
		if ( tr.Trace.visible() ) {
			tr.Trace << "have read alignment\n" << *current << std::endl;
		}
		alignments.push_back( *current );
	}
	return alignments;
} // read_general_aln

utility::vector1< SequenceAlignment > read_general_aln_file(
	std::string const & filename
) {
	utility::io::izstream input( filename.c_str() );
	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
		utility::vector1< SequenceAlignment > aligns;
		return aligns;
	} else {
		return read_general_aln( input );
	}
}

utility::vector1< SequenceAlignment > read_grishin_aln_file(
	std::string const & filename
) {
	utility::io::izstream input( filename.c_str() );
	utility::vector1< SequenceAlignment > alignments;
	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
	}

	std::string line, id1, id2, dummy;
	SequenceAlignmentOP current( new SequenceAlignment );
	while ( getline( input, line ) ) {
		if ( line.substr(0,2) == "--" ) {
			if ( current->size() > 0 ) {
				alignments.push_back( *current );
			}
			current = SequenceAlignmentOP( new SequenceAlignment );
		} else if ( line.substr(0,2) == "##" ) {
			std::istringstream line_input( line );
			line_input >> dummy >> id1 >> id2;
		} else if ( line.substr(0,1) == "#" ) {
			// do nothing
		} else if ( line.substr(0,5) == "score" ) {
			std::istringstream line_input( line );
			Real score( 0.0 );
			line_input >> dummy >> score;
			current->score( score );

			using ObjexxFCL::string_of;
			core::Size count(1);
			line_input >> score;
			while ( !line_input.fail() ) {
				count++;
				current->score( (std::string) ("score" + string_of(count)), score ) ;
				line_input >> score;
			}
		} else {
			core::Size start;
			std::string myseq, id;
			std::istringstream line_input( line );
			line_input >> start >> myseq;
			if ( ! line_input.fail() ) {
				start = start + 1; // convert from zero-based index

				if ( current->size() >= 1 ) {
					id = id2;
				} else {
					id = id1;
				}

				SequenceOP new_seq( new Sequence( myseq, id, start ) );
				current->add_sequence( new_seq );
			}
		}
	} // while( getline( input, line ) )

	if ( current->size() > 0 ) alignments.push_back( *current );

	return alignments;
}

Size
n_correctly_aligned_positions(
	SequenceAlignment & candidate_aln,
	SequenceAlignment & true_aln
) {
	runtime_assert( candidate_aln.size() == true_aln.size() );

	using core::id::SequenceMapping;
	SequenceMapping true_map      = true_aln.     sequence_mapping( 1, 2 );
	SequenceMapping candidate_map = candidate_aln.sequence_mapping( 1, 2 );

	Size n_correct(0);
	for ( Size i = 1; i <= true_aln.length(); ++i ) {
		Size resi_idx( true_aln.sequence(1)->start() + i );
		if ( true_map[ resi_idx ] == 0 ) {
			// true_map has a gap here
			if ( candidate_map[ resi_idx ] == 0 ) {
				// candidate_map also has a gap here, it's a true negative
				++n_correct;
			}
		} else {
			// true_map does not have a gap here
			if ( candidate_map[ resi_idx ] == true_map[ resi_idx ] ) {
				// candidate_map also does not have a gap, this is a true positive
				++n_correct;
			}
		}
	} // true_aln.length()

	return n_correct;
} // alignment_quality

SequenceAlignment steal_alignment(
	SequenceAlignment aln_to_steal,
	utility::vector1< SequenceOP > seqs
) {

	using std::string;
	using utility::vector1;
	using core::id::SequenceMapping;

	runtime_assert( aln_to_steal.size() == seqs.size() );

	// insertion positions
	vector1< Size > insertion_positions;
	for ( Size ii = 1; ii <= aln_to_steal.size(); ++ii ) {
		Size const start( aln_to_steal.sequence(ii)->start() );
		insertion_positions.push_back( start );
	}

	// internal gaps
	for ( Size ii = 1; ii <= aln_to_steal.length(); ++ii ) {
		for ( Size jj = 1; jj <= aln_to_steal.size(); ++jj ) {
			if ( aln_to_steal.sequence(jj)->is_gap(ii) ) {
				seqs[jj]->insert_gap( insertion_positions[jj] );
			}
			++insertion_positions[jj];
		}
	}

	// leading gaps
	for ( Size ii = 1; ii <= aln_to_steal.size(); ++ii ) {
		Size const start( aln_to_steal.sequence(ii)->start() );
		for ( Size jj = 1; jj < start; ++jj ) {
			seqs[ii]->delete_position( 1 );
		}
		seqs[ii]->start( start );
	}

	// trailing gaps
	for ( Size jj = 1; jj <= seqs.size(); ++jj ) {
		Size const seq_length( seqs[jj]->length() );
		Size const desired_length( aln_to_steal.length() );
		Size const n_to_delete( seq_length - desired_length + 1 );

		if ( seq_length < desired_length ) {
			std::cout << "--------------------------------------------------" << std::endl;
			std::cout << aln_to_steal << std::endl;
			std::cout << "seq: " << seqs[jj]->to_string() << std::endl;
			std::cout << "seqs[jj]->sequence(): " << seqs[jj]->sequence() << std::endl;
			std::cout << "length = " << seqs[jj]->length()
				<< " desired_length = " << desired_length << std::endl;
			std::cout << "n_to_delete = " << n_to_delete << std::endl;
			utility_exit_with_message( "error!" );
		}
		for ( Size ii = 1; ii < n_to_delete; ++ii ) {
			seqs[jj]->delete_position( seqs[jj]->length() );
		}
	}

	// do a quick consistency check here to make sure that sequences are
	// identical. If they're not, something has gone wrong!
	runtime_assert( aln_to_steal.size() == seqs.size() );
	for ( Size ii = 1; ii <= aln_to_steal.size(); ++ii ) {
		using std::string;
		string const aln_seq( aln_to_steal.sequence(ii)->sequence() );
		string const stolen_seq( seqs[ii]->sequence() );
		bool error( aln_seq.length() != stolen_seq.length() );

		// quick hack to deal with some sequences - X compared to
		// anything should be ok.
		if ( aln_seq != stolen_seq ) {
			for ( Size idx = 1; idx <= stolen_seq.length(); ++idx ) {
				error = ( error && aln_seq[idx] != 'X' && stolen_seq[idx] != 'X' &&
					aln_seq[idx] != stolen_seq[idx]
				);
			}
		}

		if ( error ) {
			std::string msg( "sequences are not the same!\n" );
			msg += "to_steal:\n";
			msg += aln_to_steal.sequence(ii)->to_string() + "\n";
			msg += "stolen:\n";
			msg += seqs[ii]->to_string() + "\n";
			utility_exit_with_message( msg );
		}
	}

	SequenceAlignment new_aln;
	for ( Size j = 1; j <= seqs.size(); ++j ) {
		new_aln.add_sequence( seqs[j] );
	}

	return new_aln;
} // steal_alignment

SequenceAlignment mapping_to_alignment(
	core::id::SequenceMapping const & mapping,
	SequenceOP seq1_orig,
	SequenceOP seq2_orig
) {
	SequenceOP seq1( new Sequence( "", seq1_orig->id(), seq1_orig->start() ) );
	SequenceOP seq2( new Sequence( "", seq2_orig->id(), seq2_orig->start() ) );

	runtime_assert( seq1->length() == seq1->ungapped_length() );
	runtime_assert( seq2->length() == seq2->ungapped_length() );

	core::id::SequenceMapping rev_mapping( mapping );
	rev_mapping.reverse();

	//std::cout << "forward: " << std::endl << mapping << std::endl;
	//std::cout << "reverse: " << std::endl << rev_mapping << std::endl;

	//Size idx1 = i_stop, idx2 = i_stop;
	Size ngaps1 = seq1_orig->start()-1, ngaps2 = seq2_orig->start()-1;
	Size idx1 = seq1_orig->start(), idx2 = seq2_orig->start();
	while ( idx1 <= mapping.size1() && idx2 <= rev_mapping.size1() ) {
		if ( mapping[ idx1 ] == 0 ) {
			// gap in sequence 2
			seq2->append_gap();
			seq1->append_char( (*seq1_orig)[idx1-ngaps1] );
			++idx1;
		} else if ( rev_mapping[ idx2 ] == 0 ) {
			// gap in sequence 1
			seq1->append_gap();
			seq2->append_char( (*seq2_orig)[idx2-ngaps2] );
			++idx2;
		} else {
			seq1->append_char( (*seq1_orig)[idx1-ngaps1] );
			seq2->append_char( (*seq2_orig)[idx2-ngaps2] );
			++idx1;
			++idx2;
		}
	}

	// trailing gaps
	while ( idx1-ngaps1 <= seq1_orig->length() ) {
		seq2->append_gap();
		seq1->append_char( (*seq1_orig)[idx1-ngaps1] );
		idx1++;
	}
	while ( idx2-ngaps2 <= seq2_orig->length() ) {
		seq1->append_gap();
		seq2->append_char( (*seq2_orig)[idx2-ngaps2] );
		idx2++;
	}

	runtime_assert( seq1->length() == seq2->length() );

	SequenceAlignment align;
	align.add_sequence( seq1 );
	align.add_sequence( seq2 );
	return align;
}

core::id::SequenceMapping transitive_map(
	core::id::SequenceMapping const & map1,
	core::id::SequenceMapping const & map2
) {
	//runtime_assert( map1.size2() == map2.size1() );
	//runtime_assert( map1.seq2 () == map2.seq1 () );

	// maps from sequence 1 in map1 to sequence 2 in map2
	core::id::SequenceMapping new_mapping( map1.size1(), map2.size2() );
	for ( Size ii = 1; ii <= map1.size1(); ++ii ) {
		if ( map1[ ii ] != 0 ) new_mapping[ii] = map2[ map1[ ii ] ];
	}
	return new_mapping;
}

core::id::SequenceMapping map_seq1_seq2(
	core::sequence::SequenceOP seq1,
	core::sequence::SequenceOP seq2
) {
	using namespace core::sequence;
	core::sequence::SequenceOP copy1( seq1->clone() );
	core::sequence::SequenceOP copy2( seq2->clone() );

	bool success( false ); // pessimism by default
	core::id::SequenceMapping retval;
	SWAligner sw_align;
	ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -4, -1 ) );

	// test aligning with gaps
	if ( !success ) {
		SequenceAlignment intermediate = sw_align.align( copy1, copy2, ss );
		if ( intermediate.identities() != intermediate.length() ) {
			success = false;
		} else {
			success = true;
		}
		retval = intermediate.sequence_mapping( 1, 2 );
	}

	// test aligning with no gaps
	if ( !success ) {
		copy1->sequence( seq1->ungapped_sequence() );
		copy2->sequence( seq2->ungapped_sequence() );

		SequenceAlignment intermediate = sw_align.align( copy1, copy2, ss );
		if ( intermediate.identities() != intermediate.length() ) {
			tr.Warning << "potential mismatch between sequence from alignment ";
			tr.Warning << "and sequence from PDB!" << std::endl;
			tr.Warning << "alignment: " << std::endl << intermediate
				<< std::endl;
		}
		retval = intermediate.sequence_mapping( 1, 2 );
	}

	return retval;
}

core::sequence::SequenceAlignment align_naive(
	core::sequence::SequenceOP seq1,
	core::sequence::SequenceOP seq2
) {
	using namespace core::sequence;
	core::sequence::SequenceOP copy1( seq1->clone() );
	core::sequence::SequenceOP copy2( seq2->clone() );

	SWAligner sw_align;
	ScoringSchemeOP ss( new SimpleScoringScheme( 6, 1, -2, -1 ) );

	SequenceAlignment aln = sw_align.align( copy1, copy2, ss );
	return aln;
}

core::sequence::SequenceAlignment align_poses_naive(
	core::pose::Pose & pose1,
	core::pose::Pose & pose2
) {
	using namespace core::sequence;
	SequenceOP seq1( new Sequence( pose1.sequence(), "pose1", 1 ) );
	SequenceOP seq2( new Sequence( pose2.sequence(), "pose2", 1 ) );
	return align_naive(seq1,seq2);
}


utility::vector1< Real >
get_maximum_scores(
	core::sequence::ScoringSchemeOP ss,
	core::sequence::SequenceOP seq
) {
	using core::Real;
	using utility::vector1;

	vector1< Real > scores( seq->length(), 0.0 );
	for ( Size ii = 1; ii <= seq->length(); ++ii ) {
		scores[ ii ] = ss->score( seq, seq, ii, ii );
	}

	return scores;
}

core::sequence::SequenceAlignment
alignment_from_pose(
	core::pose::Pose & pose
) {
	using core::pose::get_comment;
	using std::string;
	using core::sequence::Sequence;
	using core::sequence::SequenceOP;
	using core::sequence::SequenceAlignment;
	string q_seq, t_seq;

	bool success(
		get_comment( pose, "query_alignment   ", q_seq ) &&
		get_comment( pose, "template_alignment", t_seq )
	);

	SequenceAlignment aln;
	if ( !success ) {
		tr.Error << "Can't extract alignment from pose!" << std::endl;
		tr.Error << "query_aln:    " << q_seq << std::endl;
		tr.Error << "template_aln: " << t_seq << std::endl;
		tr.flush_all_channels();
		return aln;
	}

	SequenceOP query( new Sequence ), templ( new Sequence );
	std::istringstream q_in( q_seq ), t_in( t_seq );
	query->read_data( q_in );
	templ->read_data( t_in );

	aln.add_sequence( query );
	aln.add_sequence( templ );
	tr.Debug << "extracted sequence alignment from Pose: " << std::endl
		<< aln << std::endl;
	return aln;
}

void alignment_into_pose(
	core::sequence::SequenceAlignment const & aln,
	core::pose::Pose & pose
) {
	add_comment( pose, "query_alignment   ", aln.sequence(1)->to_string() );
	add_comment( pose, "template_alignment", aln.sequence(2)->to_string() );
}

core::Real
calpha_superimpose_with_mapping(
	core::pose::Pose & mod_pose,
	core::pose::Pose const & ref_pose,
	core::id::SequenceMapping const & mapping // mod_pose -> ref_pose
) {
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::AtomID::BOGUS_ATOM_ID() );
	Size const mod_resn( mod_pose.size() );
	Size const ref_resn( ref_pose.size() );
	static std::string const atom_name("CA");
	for ( Size mod_resi = 1; mod_resi <= mod_resn; ++mod_resi ) {
		Size const ref_resi( mapping[mod_resi] );
		if ( ref_resi && mod_resi <= mod_resn && ref_resi <= ref_resn ) {
			if ( ! mod_pose.residue(mod_resi).has(atom_name) ) continue;
			if ( ! ref_pose.residue(ref_resi).has(atom_name) ) continue;

			id::AtomID const id1( mod_pose.residue(mod_resi).atom_index(atom_name), mod_resi );
			id::AtomID const id2( ref_pose.residue(ref_resi).atom_index(atom_name), ref_resi );
			atom_map.set( id1, id2 );
		}
	}
	return core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}

////////////////////////////////////////////////////////
utility::vector1< Size >
strip_spacers( std::string & sequence, bool const annotations_in_brackets /* = true */ )
{
	std::string new_sequence;
	utility::vector1< Size > spacer_pos;
	Size offset = 0;
	for ( Size n = 1; n <= sequence.size(); n++ ) {
		if ( annotations_in_brackets ) {
			if ( sequence[ n - 1 ] == '[' ) {
				new_sequence += sequence[ n - 1 ];
				++n;
				++offset;
				while ( sequence[ n - 1 ] != ']' ) {
					new_sequence += sequence[ n - 1 ];
					++n;
					++offset;
				}
			}
			if ( sequence[ n - 1 ] == '[' ) {
				new_sequence += sequence[ n - 1 ];
			}
		}
		if ( spacers.has_value( sequence[ n - 1 ] ) ) {
			Size spacer_position = new_sequence.size() - offset;
			if ( new_sequence.size() > 0 ) spacer_pos.push_back( spacer_position );
			continue;
		}
		new_sequence += sequence[ n - 1 ];
	}
	sequence = new_sequence;
	return spacer_pos;
}


///////////////////////////////////////////////////////////////////////////////////////
// Converts sequences to one-letter-sequences, but outputs a map of any fullnames (as would be enclosed in brackets, like for Z[IGU]).
std::map< Size, std::string >
parse_out_non_standard_residues( vector1< core::sequence::SequenceOP > & fasta_sequences ) {

	using namespace core::sequence;
	std::map< Size, std::string > non_standard_residues;
	vector1< core::sequence::SequenceOP > fasta_sequences_new;

	Size offset( 0 );
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		std::string sequence = fasta_sequences[ n ]->sequence();
		std::map< Size, std::string > non_standard_residues_local = parse_out_non_standard_residues( sequence );

		SequenceOP new_sequence( new Sequence( sequence, fasta_sequences[n]->id()) );
		new_sequence->spacer_positions( fasta_sequences[n]->spacer_positions() );
		fasta_sequences_new.push_back( new_sequence );
		for ( auto & it : non_standard_residues_local ) {
			non_standard_residues[ it.first + offset ] = it.second;
		}

		offset += sequence.size();
	}

	fasta_sequences = fasta_sequences_new;
	return non_standard_residues;
}

///////////////////////////////////////////////////////////////////////////////////////
std::map< Size, std::string >
parse_out_non_standard_residues( std::string & sequence ) {
	utility::vector1< std::string > fullname_list; // a vector of non-standard full names
	std::vector< Size > oneletter_to_fullname_index; // for each one-letter sequence, zero means no fullname given
	std::string one_letter_sequence;

	core::pose::parse_sequence( sequence, fullname_list, oneletter_to_fullname_index, one_letter_sequence );
	sequence = one_letter_sequence;

	std::map< Size, std::string > non_standard_residues;
	for ( Size k = 1; k <= oneletter_to_fullname_index.size(); k++ ) {
		Size const pos = oneletter_to_fullname_index[ k - 1 ];
		if ( pos > 0 ) non_standard_residues[ k ] = fullname_list[ pos ];
	}
	return non_standard_residues;
}

/// @brief Convert sequence string to fasta string with only 80
/// characters per line
std::string
convert_to_fasta( std::string const & pname, std::string const & seq )
{
	std::string fasta_str( ">" + pname + "\n" );
	for ( core::Size i=0; i<seq.size(); i++ ) {
		fasta_str += seq[i];
		if ( (i+1) % 80 == 0 ) {
			fasta_str += "\n";
		} //if 80 chars
	} // for each residue
	return fasta_str;
}

/// @brief Create fasta file from sequence string. Differs from
/// output_fasta_file in that the output is in valid FASTA format,
/// rather than including NTerm/CTerm tags.
std::string
create_fasta_file( std::string const & pname, std::string const & seq )
{
	std::string fasta_filename( pname + ".fasta" );
	// get rid of the path and just use the current directory
	if ( fasta_filename.find('/') != std::string::npos ) {
		fasta_filename = fasta_filename.substr( fasta_filename.find_last_of( '/' )+1, std::string::npos );
	}
	std::string fasta( convert_to_fasta( pname, seq ) );
	tr.Debug << "Fasta: " << fasta << std::endl;
	tr << "Fasta filename: " << fasta_filename << std::endl;

	// Write out a fasta file
	std::ofstream fastafile( fasta_filename.c_str() );
	if ( !fastafile.is_open() ) {
		utility_exit_with_message( "Could not open " + fasta_filename + " for writing." );
	}
	fastafile << fasta << std::endl;
	if ( fastafile.bad() ) {
		utility_exit_with_message( "Encountered an error writing to " + fasta_filename );
	}
	fastafile.close();
	return fasta_filename;
}

// @brief output annotated_sequence (as well as numbering/chain info) to FASTA-like file.
void
output_fasta_file( std::string const & fasta_filename, core::pose::Pose const & pose  )
{
	// Write out a fasta file
	std::ofstream fastafile( fasta_filename.c_str() );
	if ( !fastafile.is_open() ) {
		utility_exit_with_message( "Error: could not open " + fasta_filename + " for writing." );
	}

	utility::vector1< int > res_vector;
	utility::vector1< char > chain_vector;
	utility::vector1< std::string > segid_vector;
	for ( Size i = 1; i <= pose.size(); i++ ) {
		res_vector.push_back(   pose.pdb_info()->number( i ) );
		chain_vector.push_back( pose.pdb_info()->chain( i ) );
		segid_vector.push_back( pose.pdb_info()->segmentID( i ) );
	}
	fastafile << "> " << tag_from_pose( pose ) << " " << make_tag_with_dashes( res_vector, chain_vector, segid_vector ) << std::endl;

	fastafile << pose.annotated_sequence() << std::endl;

	if ( fastafile.bad() ) {
		utility_exit_with_message( "Encountered an error writing to " + fasta_filename );
	}
	fastafile.close();
}

} // sequence
} // core
