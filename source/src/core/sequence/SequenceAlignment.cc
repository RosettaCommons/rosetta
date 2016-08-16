// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SequenceAlignment.hh
/// @brief class definition for a multiple sequence alignment
/// @author James Thompson

#include <core/types.hh>

#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>

#include <ObjexxFCL/string.functions.hh>

#include <iostream>
#include <string>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace sequence {

/// @details moving from .cc file
SequenceAlignment::~SequenceAlignment() {
	// clear(); // APL NOTE: originally, this dstor called clear(), but clear() doesn't do anything that the dstor itself doesn't do
}

static THREAD_LOCAL basic::Tracer tr( "core.sequence.SequenceAlignment" );

void SequenceAlignment::add_sequence( SequenceOP myseq ) {
	sequences_.push_back( myseq );
}

Size SequenceAlignment::size() const {
	return sequences_.size();
}

Size SequenceAlignment::length() const {
	if ( size() < 1 ) return 0;
	else return sequences_[1]->length();
}

SequenceOP SequenceAlignment::sequence( Size idx ) const {
	if ( idx > size() ) {
		using ObjexxFCL::string_of;
		std::string msg("");
		msg += "Requested sequence " + string_of(idx)
			+ " but alignment only has " + string_of(size()) + " sequences\n";
		msg += "Alignments:\n" + to_string() + "\n";
		utility_exit_with_message( msg );
		//runtime_assert( idx <= size() );
	}
	return sequences_[idx];
}

std::map< std::string, core::Real > SequenceAlignment::scores() const {
	return scores_;
}

void SequenceAlignment::scores( std::map< std::string, core::Real >  new_scores ) {
	scores_ = new_scores;
}

Real SequenceAlignment::score() const {
	static std::string const name("score");
	return score(name);
}

void SequenceAlignment::score( Real const & sc ) {
	score( "score", sc );
}

Real SequenceAlignment::score( std::string const & name ) const {

	std::map< std::string, core::Real >::const_iterator it(
		scores_.find(name)
	);

	if ( it == scores_.end() ) return 0;
	return it->second;
}

void SequenceAlignment::score( std::string const & name, Real const value ) {
	scores_[name] = value;
}


std::string SequenceAlignment::to_string() const {
	std::string retval("");
	retval += "# score " + ObjexxFCL::string_of(score()) + "\n";
	// add comments here
	for ( Size i = 1; i <= size(); ++i ) {
		retval += sequence(i)->to_string() + "\n";
	}
	retval += "--\n";
	return retval;
}

std::string SequenceAlignment::alignment_id() const {
	return ObjexxFCL::uppercased( sequence(2)->id() );
}

/// @brief super-general alignment format reader. This will never ever ever
/// have to be improved or extended.
void SequenceAlignment::read_from_file( std::string const & filename ) {
	utility::io::izstream data( filename );
	if ( !data ) {
		utility_exit_with_message(
			"ERROR: Unable to open alignment file: " +  filename
		);
	}

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream input( line );
		SequenceOP temp_seq( new Sequence );
		temp_seq->read_data( input );
		add_sequence( temp_seq );
	}
} // read_from_file

void SequenceAlignment::read_data( std::istream & in ) {
	std::string line;
	while ( getline( in, line ) ) {
		std::istringstream input( line );
		if ( line.substr(0,1) == "#" ) {
			if ( line.substr(2,5) == "score" ) {
				std::string dummy;
				Real sc;
				input >> dummy >> dummy >> dummy >> sc;
				score( sc );
			}
		} else if ( line.substr(0,2) == "--" ) {
			break;
		}

		SequenceOP temp_seq( new Sequence );
		temp_seq->read_data( input );
		add_sequence( temp_seq );
	}
} // read_data

/// @brief returns a SequenceMapping of the sequence at index idx1 mapped to
/// the sequence at idx2.
core::id::SequenceMapping SequenceAlignment::sequence_mapping(
	Size const idx1,
	Size const idx2
) const {
	runtime_assert( idx1 <= size() && idx1 > 0 );
	runtime_assert( idx2 <= size() && idx2 > 0);

	using core::id::SequenceMapping;

	utility::vector1< SequenceOP >::const_iterator
		it1( sequences_.begin() + (idx1 - 1) ),
		it2( sequences_.begin() + (idx2 - 1) );

	SequenceMapping mapping( (*it1)->ungapped_length(), (*it2)->ungapped_length() );

	for ( Size pos = 1; pos <= length(); ++pos ) {
		Size const seq1_pos( (*it1)->resnum(pos) );
		Size const seq2_pos( (*it2)->resnum(pos) );
		if ( seq1_pos != 0 && seq2_pos != 0 ) {
			mapping[ seq1_pos ] = seq2_pos;

			//char const aa1( (*it1)->sequence().at(seq1_pos-1) );
			//char const aa2( (*it2)->sequence().at(seq2_pos-2) );
			//std::cout << "texdebug: adding mapping of "
			// << seq1_pos << "," << aa1 << " => "
			// << seq2_pos << "," << aa2 << std::endl;
		}
	}

	return mapping;
} // sequence_mapping

void SequenceAlignment::remove_gapped_positions() {
	Size pos(1);
	while ( pos <= length() ) {
		// delete this new position if the entire column is gapped.
		bool delete_column( true );
		for ( utility::vector1< SequenceOP >::iterator it = sequences_.begin(),
				end = sequences_.end();
				it != end; ++it
				) {
			if ( !(*it)->is_gap(pos) ) {
				// std::cout << "not deleting column because sequence " << *it
				// << " has no gap at position " << pos << std::endl;
				delete_column = false;
			}
		} // for sequences

		if ( delete_column ) {
			for ( utility::vector1< SequenceOP >::iterator it = sequences_.begin(), end = sequences_.end();
					it != end; ++it
					) {
				// std::cout << "deleting column " << pos << " from " << *it << std::endl;
				(*it)->delete_position( pos );
			}
		} else {
			++pos;
		}

	} // while pos <= length()
} // remove_gapped_positions

Real SequenceAlignment::calculate_score_sum_of_pairs(
	ScoringSchemeOP ss
) const {
	using core::Real;
	using utility::vector1;

	Real score( 0.0 );
	vector1< Real > scores = calculate_per_position_scores( ss );
	// something mis-understood about std::accumulate?
	//std::accumulate( scores.begin(), scores.end(), 0 );
	typedef vector1< Real >::const_iterator iter;
	for ( iter it = scores.begin(), end = scores.end(); it != end; ++it ) {
		score += *it;
	}

	return score;
} // calculate_score_sum_of_pairs

utility::vector1< Real > SequenceAlignment::calculate_per_position_scores(
	ScoringSchemeOP ss
) const {
	utility::vector1< Real > scores( length(), 0.0 );
	for ( Size i = 1; i <= size(); ++i ) {
		for ( Size j = i + 1; j <= size(); ++j ) {
			bool last_pos_was_gapped( false );
			for ( Size k = 1; k <= length(); ++k ) {
				if ( sequence(i)->is_gap(k) && sequence(j)->is_gap(k) ) {
					// maybe rethink this? you could argue for a positive score as well,
					// such as 2 * gap_open or 2 * gap_extend, depending on
					// last_pos_was_gapped.
					last_pos_was_gapped = true;
					scores[k] += 0.0;
				} else if ( sequence(i)->is_gap(k) || sequence(j)->is_gap(k) ) {
					if ( last_pos_was_gapped ) scores[k] += ss->gap_extend();
					else                       scores[k] += ss->gap_open();
					last_pos_was_gapped = true;
				} else {
					bool over_run(false);
					if ( sequences_[i]->length() < k ) {
						over_run = true;
						tr.Error << "Attempting to run off the end of sequence!" << std::endl;
						tr.Error << "asked for position " << k << " in seq:" << std::endl;
						tr.Error << sequences_[i]->to_string() << std::endl;
					}
					if ( sequences_[j]->length() < k ) {
						over_run = true;
						tr.Error << "Attempting to run off the end of sequence!" << std::endl;
						tr.Error << "asked for position " << k << " in seq:" << std::endl;
						tr.Error << sequences_[j]->to_string() << std::endl;
					}
					if ( !over_run ) {
						scores[k] += ss->score( sequences_[i], sequences_[j], k, k );
						last_pos_was_gapped = false;
					}
				}
			} // col k
		} // seq j
	} // seq i

	return scores;
}

void SequenceAlignment::data_integrity_check() const {
	for ( Size jj = 1; jj <= size(); ++jj ) {
		if ( sequences_[jj]->length() != length() ) {
			std::string msg( "Error: length mismatch between sequence and alignment" );
			msg += "problem with sequence: " + sequences_[jj]->to_string();
			msg += "alignment: " + to_string();
			utility_exit_with_message( msg );
		}
	} // sequence jj
}

Size SequenceAlignment::identities() const {
	Size n_ident(0);

	for ( Size i = 1; i <= length(); ++i ) {
		bool ident(true);
		data_integrity_check();

		for ( Size j = 2; j <= sequences_.size(); ++j ) {
			if ( (*sequences_[j])[i] != (*sequences_[1])[i] ) {
				ident = false;
			}
		} // for sequences

		if ( ident ) ++n_ident;
	} // for ( size i = 1 )
	return n_ident;
} // identities

Size SequenceAlignment::gapped_positions() const {
	Size n_gapped(0);

	for ( Size ii = 1; ii <= length(); ++ii ) {
		if ( is_gapped(ii) ) ++n_gapped;
	} // for ( size i = 1 )
	return n_gapped;
} // gapped_positions

bool SequenceAlignment::sequences_are_in_frame() const {
	bool in_frame( true );

	for ( Size i = 1; i <= length(); ++i ) {
		utility::vector1< SequenceOP >::const_iterator seq1, it, end;

		seq1 = sequences_.begin();
		if ( (*seq1)->is_gap( i ) ) continue;

		for ( it = seq1 + 1, end = sequences_.end(); it != end; ++it ) {
			if ( (*it)->is_gap( i ) ) continue;

			if ( (*seq1)->resnum(i) != (*it)->resnum(i) ) {
				in_frame = false;
				break; // out of inner loop
			}
		}

		if ( in_frame == false ) break; // out of outer loop
	}

	return in_frame;
}

utility::vector1< core::Size > SequenceAlignment::sequence_indices(
	core::Size const column
) const {
	runtime_assert( column <= length() );
	runtime_assert( column > 0 );

	utility::vector1< core::Size > indices;

	typedef utility::vector1< SequenceOP > seqlist;
	for ( seqlist::const_iterator it = sequences_.begin(), end = sequences_.end();
			it != end; ++it
			) {
		indices.push_back( (*it)->resnum(column) );
	}

	return indices;
}


void SequenceAlignment::comment( std::string const & comment ) {
	comments_.push_back( comment );
}

utility::vector1< std::string > SequenceAlignment::comments() const {
	return comments_;
}

Real SequenceAlignment::max_gap_percentage() const {
	typedef utility::vector1< SequenceOP > seqlist;
	Real max_gp(0.0);
	for ( seqlist::const_iterator it = sequences_.begin(), end = sequences_.end();
			it != end; ++it
			) {
		Real gap_percentage = static_cast< Real >
			( (*it)->length() - (*it)->ungapped_length() );
		gap_percentage = gap_percentage / static_cast< Real > ( length() );
		max_gp = std::max( max_gp, gap_percentage );
		//std::cout << "length = " << (*it)->length() << std::endl;
		//std::cout << "ungapped_length = " << (*it)->ungapped_length() << std::endl;
	}

	return max_gp;
}

bool SequenceAlignment::is_gapped( Size const col_idx ) const {
	bool gapped( false );
	for ( Size jj = 1; jj <= size(); ++jj ) {
		if ( sequences_[jj]->is_gap( col_idx ) ) {
			gapped = true;
			break;
		}
	} // for sequences
	return gapped;
}

//void SequenceAlignment::trim_terminal_gaps() {
//  // find the first non-gap character
// utility::vector1< core::Size > nterm_gaps;
// typedef utility::vector1< SequenceOP > seqlist;
// for ( seqlist::const_iterator it = sequences_.begin(), end = sequences_.end();
//  it != end; ++it
// ) {
//
// }
//} // trim_terminal_gaps

std::ostream & operator << (
	std::ostream & out,
	const SequenceAlignment & sa
) {
	out << "score: " << sa.score()
		<< " identities: " << sa.identities()
		<< "/" << sa.length()
		<< " gaps: " << sa.gapped_positions()
		<< "/" << sa.length()
		<< std::endl;

	for ( Size i = 1; i <= sa.size(); ++i ) {
		out << (*sa.sequence(i)) << std::endl;
	}
	return out;
}


std::istream & operator>> (
	std::istream & in,
	SequenceAlignment & aln
) {
	aln.read_data( in );
	return in;
}

SequenceAlignment::SequenceAlignment(
	SequenceAlignment const & src
) :
	ReferenceCount ( src )
{
	*this = src;
}

SequenceAlignment &
SequenceAlignment::operator =(
	SequenceAlignment const & src
) {
	clear();

	scores_ = src.scores();

	// copy sequences, be sure to call clone manually
	for ( Size i = 1; i <= src.size(); ++i ) {
		core::sequence::SequenceOP copy
			= src.sequence(i)->clone();
		add_sequence( copy );
	}

	comments_ = src.comments();
	return *this;
}

bool
operator<(
	SequenceAlignment const & lhs, SequenceAlignment const & rhs
) {
	// compare number of sequences first
	if ( lhs.size() < rhs.size() ) { return true; }

	// if sequences are the same, compare concatentation of
	// all sequences from each alignment lexicographically
	std::string l_str(""), r_str("");
	for ( Size ii = 1; ii <= lhs.size(); ++ii ) {
		l_str += lhs.sequence(ii)->sequence();
		r_str += rhs.sequence(ii)->sequence();
	}

	return ( l_str < r_str );
}

void SequenceAlignment::printGrishinFormat (
	std::ostream & out
) const
{

	//if ( size() < 2 ) return;

	out << "## " << sequence(1)->id() << " " << sequence(2)->id() << std::endl;
	out << "#  " << std::endl;
	//out << "scores_from_program: 0.000000 " << score() << std::endl;
	out << "scores_from_program:";
	// print scores in order of sorted keys
	using std::map;
	using std::string;
	using utility::vector1;

	vector1< string > keys;
	for ( map< string, Real >::const_iterator it = scores_.begin(), end = scores_.end(); it != end; ++it ) {
		keys.push_back( it->first );
	}
	std::sort( keys.begin(), keys.end() );

	for ( vector1< string >::const_iterator it = keys.begin(), end = keys.end(); it != end; ++it ) {
		out << " " << score(*it);
		//std::cout << "score(" << *it << "," << score(*it) << ")" << std::endl;
	}
	out << std::endl;

	// print out sequences
	for ( Size i = 1; i <= size(); ++i ) {
		out << sequence(i)->start()-1 << " " << (*sequence(i)).sequence() << std::endl;
	}
	out << "--" << std::endl;
}

} // sequence
} // core
