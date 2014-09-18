// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/sequence/util.cc
/// @brief small bundle of utilities for dealing with sequences
/// @author James Thompson

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
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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

static thread_local basic::Tracer tr( "core.sequence" );

void read_all_alignments(const std::string& format,
                         const utility::vector1<std::string>& files,
                         utility::vector1<SequenceAlignment>* alignments) {
  using std::string;
  using utility::vector1;

  assert(alignments);
  for (vector1<string>::const_iterator i = files.begin(); i != files.end(); ++i) {
    vector1<SequenceAlignment> current = read_aln(format, *i);
    std::copy(current.begin(), current.end(), std::back_inserter(*alignments));
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

	basic::T("core.sequence.DerivedSequenceMapping") << "align1: " << align1 << "\nalign2: " << align2 <<
		"\nseq1: " << seq1 << "\nseq2: " << seq2 << '\n';

	runtime_assert( mapping.size1() == seq1.size() );
	mapping.size2( seq2.size() ); // set sequence 2 size
} // read_alignment_file

vector1< string > read_fasta_file_str( std::string const & filename ) {
	vector1< SequenceOP > seqs = read_fasta_file( filename );
	vector1< string > seq_strings;
	for ( vector1< SequenceOP >::const_iterator it = seqs.begin(), end = seqs.end(); it != end; ++it ) {
		seq_strings.push_back( (*it)->sequence() );
	}
	return seq_strings;
}

vector1< SequenceOP > read_fasta_file( std::string const & filename ) {
	vector1< SequenceOP > sequences;

	utility::io::izstream input( filename.c_str() );
	std::string line, current_sequence = "", current_id = "empty";

	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
		return sequences;
	}

	while( getline( input, line ) ) {
		if ( line.substr(0,1) == ">" ) {
			if ( current_sequence != "" ) {
				ObjexxFCL::strip_whitespace( current_id );
				ObjexxFCL::strip_whitespace( current_sequence );
				sequences.push_back( new Sequence( current_sequence, current_id ) );
				current_sequence = "";
			}
			current_id = line.substr(1,line.size());
			continue;
		}
		current_sequence = current_sequence + ObjexxFCL::rstrip(line);
	}
	if ( current_sequence != "" ) {
		ObjexxFCL::strip_whitespace( current_id );
		ObjexxFCL::strip_whitespace( current_sequence );
		sequences.push_back( new Sequence( current_sequence, current_id ) );
	}

	return sequences;
} // read_fasta_file


std::string read_fasta_file_return_str( std::string const & filename ) {
	utility_exit_with_message(
		"This function is redundant with the functions above it. Ask for help with C++ if you need it, but this is really embarassing duplication."
	);
	std::string sequences;

	utility::io::izstream input( filename.c_str() );

	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
		return sequences;
	}

	std::string line;
	char aa;

	// the following line is dangerous, because FASTA files don't necessarily
	// have a > line.
	getline( input, line ); // skip the > line
	while( getline( input,line) ) {
		std::istringstream line_stream( line );
		while( line_stream >> aa ) {
			sequences += aa;
		}
	}

	return sequences;
} // read_fasta_file


core::sequence::DerivedSequenceMapping simple_mapping_from_file( std::string const & filename ) {
	// file I/O stuff
	utility::io::izstream input( filename.c_str() );
	std::string line, current_sequence = "";

	if ( !input ) {
		utility_exit_with_message( "Warning: can't open file " + filename + "!" );
	}

	utility::vector1< std::pair< Size,Size > > aligned;
	Size max_resi = 0, max_resj = 0;

	std::string seq1, seq2;
	Size start_seq2( 0 );
	while( getline( input, line ) ) {
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

	//	if ( tr.Trace.visible() ) mapping.show( tr.Trace );

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
		typedef vector1< FileName >::const_iterator iter;
		for ( iter it = fns.begin(), end = fns.end(); it != end; ++it ) {
			vector1< SequenceOP > temp_seqs( read_fasta_file( *it ) );
			for ( vector1< SequenceOP >::const_iterator s_it = temp_seqs.begin(),
						s_end = temp_seqs.end(); s_it != s_end; ++s_it
			) {
				seqs.push_back( *s_it );
			}
		}
	}

	if ( option[ in::file::pssm ].user() ) {
		vector1< FileName > fns( option[ in::file::pssm ]() );
		for ( vector1< FileName >::const_iterator it = fns.begin(), end = fns.end();
					it != end; ++it
		) {
			SequenceProfileOP prof( new SequenceProfile );
			prof->read_from_file( *it );
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
		retval = read_general_aln_file(	filename );
	} else if ( format == "grishin" ) {
		retval = read_grishin_aln_file(	filename );
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
	while( getline( input, line ) ) {
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
			current = new SequenceAlignment;
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

	std::string line, id1, id2, method, dummy;
	SequenceAlignmentOP current( new SequenceAlignment );
	while( getline( input, line ) ) {
		if ( line.substr(0,2) == "--" ) {
			if ( current->size() > 0 ) {
				alignments.push_back( *current );
			}
			current = new SequenceAlignment;
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
				std::string debug_s( seqs[jj]->sequence() );
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
			tr.Warning << "Error: potential mismatch between sequence from alignment ";
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
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
	Size const mod_resn( mod_pose.total_residue() );
	Size const ref_resn( ref_pose.total_residue() );
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

} // sequence
} // core
