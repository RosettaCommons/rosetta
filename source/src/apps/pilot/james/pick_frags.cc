// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file pick_frags.cc
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>

#include <basic/options/option.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/L1ScoringScheme.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/ScoringSchemeFactory.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
//using namespace basic;
//using namespace core::sequence;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

class VallReader {

public:
	VallReader( utility::file::FileName const & fn ) {
		read_lines( fn );
	} // ctor

	std::string get_line( Size idx ) const {
		assert( lines_.size() >= idx );
		return lines_[ idx ];
	}

	Size size() const {
		return lines_.size();
	}

	core::sequence::SequenceOP make_nmer_profile( Size idx, Size size, bool & err ) {
		std::string aa_seq;
		utility::vector1< utility::vector1< core::Real > > profile;

		// order of amino acid probabilities in the vall file
		static utility::vector1< core::chemical::AA > order;
		order.resize( 20 );
		order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
		order[ 2] = core::chemical::aa_from_oneletter_code( 'C' );
		order[ 3] = core::chemical::aa_from_oneletter_code( 'D' );
		order[ 4] = core::chemical::aa_from_oneletter_code( 'E' );
		order[ 5] = core::chemical::aa_from_oneletter_code( 'F' );
		order[ 6] = core::chemical::aa_from_oneletter_code( 'G' );
		order[ 7] = core::chemical::aa_from_oneletter_code( 'H' );
		order[ 8] = core::chemical::aa_from_oneletter_code( 'I' );
		order[ 9] = core::chemical::aa_from_oneletter_code( 'K' );
		order[10] = core::chemical::aa_from_oneletter_code( 'L' );
		order[11] = core::chemical::aa_from_oneletter_code( 'M' );
		order[12] = core::chemical::aa_from_oneletter_code( 'N' );
		order[13] = core::chemical::aa_from_oneletter_code( 'P' );
		order[14] = core::chemical::aa_from_oneletter_code( 'Q' );
		order[15] = core::chemical::aa_from_oneletter_code( 'R' );
		order[16] = core::chemical::aa_from_oneletter_code( 'S' );
		order[17] = core::chemical::aa_from_oneletter_code( 'T' );
		order[18] = core::chemical::aa_from_oneletter_code( 'V' );
		order[19] = core::chemical::aa_from_oneletter_code( 'W' );
		order[20] = core::chemical::aa_from_oneletter_code( 'Y' );

		err = false;
		Size last_resi( 1 );
		std::string last_id( "" );
		for ( Size i = idx; i <= idx + size - 1; ++i ) {
			std::string line = get_line( i );
			std::istringstream ls( line );

			Size resi;
			char ss;
			std::string id, aa, dummy;
			Real x, y, z, phi, psi, omega;
			//			utility::vector1< Real > profile_row( 0.0, order.size() );
			utility::vector1< Real > profile_row;
			profile_row.resize( order.size() );
			ls 	>> id >> aa >> ss >> resi
					>> dummy >> dummy
					>> x >> y >> z
					>> phi >> psi >> omega
					>> dummy >> dummy >> dummy >> dummy;
			aa_seq += aa;

			Size index(1);
			Real aa_prob;
			ls >> aa_prob;
			while ( !ls.fail() ) {
				//profile_row.push_back( aa_prob );
				profile_row[ order[index] ] = aa_prob;
				ls >> aa_prob;
				++index;
			}

			if ( i == idx ) { // hack for error checking
				last_id   = id;
				last_resi = resi - 1;
			}

			// error-checking.
			if ( last_resi != resi - 1 ) {
				//std::cout << "error 1: " << last_resi << " != " << resi - 1 << std::endl;
				err = true;
				break;
			}
			if ( last_id != id ) {
				//std::cout << "error 2: " << last_id << " != " << id << std::endl;
				err = true;
				break;
			}
			//if ( err ) break;

			profile.push_back( profile_row );

			// error checks for next iteration through loop
			last_resi = resi;
			last_id   = id;
			//std::cout << "profile_row.size() = " << profile_row.size() << std::endl;
		} // for ( Size i )

		core::sequence::SequenceOP seq_profile(
			new core::sequence::SequenceProfile( profile, aa_seq, "empty" )
		);
		seq_profile->sequence( aa_seq );

		return seq_profile;
	} // make_nmer_profile

	utility::vector1< std::string > make_fragment_lines(
		Size idx,
		Size size
	) {

		utility::vector1< std::string > lines;

		for ( Size i = idx; i <= idx + size - 1; ++i ) {
			std::string line = get_line( i );
			std::istringstream ls( line );

			Size resi;
			char ss;
			std::string id, aa, dummy;
			Real x, y, z, phi, psi, omega;
			utility::vector1< Real > profile_row;
			ls 	>> id >> aa >> ss >> resi
					>> dummy >> dummy
					>> x >> y >> z
					>> phi >> psi >> omega
					>> dummy >> dummy >> dummy >> dummy;

			Real aa_prob;
			ls >> aa_prob;
			while ( !ls.fail() ) {
				profile_row.push_back( aa_prob );
				ls >> aa_prob;
			}
			std::string pdbid = id.substr(0,4);
			std::string chain = id.substr(4,1);

			std::string frag_line =
				+ " " + pdbid + " " + chain + " "
				+ I( 5, resi ) + " "
				+ aa + " " + ss
				+ F( 9, 3, phi   )
				+ F( 9, 3, psi   )
				+ F( 9, 3, omega );
			lines.push_back( frag_line );
		} // for ( Size i = idx

		return lines;
	} // make_fragment_lines

private:
	void read_lines( FileName const & fn ) {
		lines_.clear();
    utility::io::izstream input( fn );
		if ( !input ) {
			utility_exit_with_message( "ERROR: Unable to open file (" + fn.name() + ")" );
		}

		std::string line;
		while ( getline( input, line ) ) {
			lines_.push_back( line );
		}
	} // get_lines

	utility::vector1< std::string > lines_;
}; // VallReader

class VallIndex {
public:
	VallIndex( Size idx, Real value ) : index_( idx ), value_( value ) {}

	Size index() const {
		return index_;
	}

	Real value() const {
		return value_;
	}

	VallIndex & operator = ( VallIndex const & rhs ) {
		if ( this == &rhs ) return *this;

		value_ = rhs.value();
		index_ = rhs.index();
		return *this;
	}

	friend bool operator < ( VallIndex const & lhs, VallIndex const & rhs ) {
		return ( lhs.value() < rhs.value() );
	}

private:
	Size index_;
	Real value_;
}; // VallIndex


struct compare_VallIndex_lt {
  bool operator()( VallIndex const & a, VallIndex const & b ) {
		return( a.value() < b.value() );
  }
};

struct compare_VallIndex_gt {
  bool operator()( VallIndex const & a, VallIndex const & b ) {
		return( a.value() > b.value() );
  }
};

// Silly little class to hold a sorted list of values. Handles the logic of
// only holding onto maxN values using a custom comparator.
template < typename T, class Comparator >
class HeapContainer {
public:
	/// @brief default constructor for HeapContainer. Sets max_values to 100.
	HeapContainer() : max_values_( 100 ), is_sorted_( false ) {}
	/// @brief constructor.
	HeapContainer( Size const max ) : max_values_( max ), is_sorted_( false ) {}

	/// @brief Returns the maximum number of values that this HeapContainer will
	/// hold.
	Size max_values() const {
		return max_values_;
	}

	/// @brief Sets the maximum number of values that this HeapContainer will
	/// hold.
	void max_values( Size new_max ) {
		max_values_ = new_max;
	}

	/// @brief Add a value to his HeapContainer. If it's greater than any of the
	/// values already inside, it will be added to the HeapContainer. Only
	/// max_values() values are stored.
	void add_value( T val ) {
		if ( values_.size() + 1 <= max_values_ ) {
			// easy case, we're under the Size limit. no need to sort.
			values_.push_back( val );
			is_sorted_ = false;
			return;
		}

		// if we've gotten here, we're full up to max_values_
		assert( values_.size() == max_values() );

		sort_myself(); // make certain that I'm sorted, set worst_value

		// no need to do anything, new value isn't large enough
		//if ( val.value() < worst_value() ) return;
		if ( Comparator()( worst(), val ) ) return;

		values_[ values_.size() ] = val;
		is_sorted_ = false;
	} // add_value

	/// @brief Returns the number of values currently in this object.
	Size size() const {
		return values_.size();
	}

	/// @brief Access the value at index idx.
	T operator [] ( Size idx ) {
		assert( idx <= values_.size() );
		return values_[idx];
	}

private:
	inline T worst() {
		sort_myself();
		return values_[ values_.size() ];
	}

	inline void sort_myself() {
		if ( !is_sorted_ ) {
			std::sort( values_.begin(), values_.end(), Comparator() );
			is_sorted_ = true;
		}
	}

	core::Size max_values_;
	bool is_sorted_;
	utility::vector1< T > values_;
}; // HeapContainer

int
main( int argc, char* argv [] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// options, random initialization
	devel::init( argc, argv );

	std::string name( option[ in::file::vall ]() );
	FileName fn( name );
	VallReader vall( name );

	// query profile
	core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );
	q_prof->read_from_checkpoint( option[ in::file::fasta ]()[1] );

	// scoring scheme for picking fragments
	std::string const scoring_scheme_type( option[ frags::scoring::profile_score ]() );
	core::sequence::ScoringSchemeFactory ssf;
	core::sequence::ScoringSchemeOP ss( ssf.get_scoring_scheme( scoring_scheme_type ) );

	std::cout << "picking fragments for query sequence: "
						<< q_prof->sequence() << std::endl;
	Size const frags_per_position( option[ frags::n_frags ]() );

	utility::vector1< Size > frag_sizes;
	frag_sizes.push_back( 1 );
	frag_sizes.push_back( 3 );
	frag_sizes.push_back( 9 );

	utility::vector1< core::Real > dummy2( vall.size(), -99.0 );
	utility::vector1< utility::vector1< core::Real > > single_per_pos_scores( q_prof->length(), dummy2 );

	typedef utility::vector1< Size > sizelist;
	for ( sizelist::const_iterator it = frag_sizes.begin(), end = frag_sizes.end();
				it != end; ++it
	) {

		Size const frag_size( *it );
		HeapContainer< VallIndex, compare_VallIndex_lt > dummy( frags_per_position );
		utility::vector1< HeapContainer< VallIndex, compare_VallIndex_lt > >
			per_position_scores( q_prof->length() - frag_size + 1, dummy );

		std::cerr << "picking fragments of length " << frag_size
							<< " in vall of size " << vall.size() << "." << std::endl;


		bool vall_err( false );
		PROF_START( basic::SEQUENCE_COMPARISON );

		for ( Size i = 1, vall_size = vall.size(); i <= vall_size - frag_size + 1; ++i ) {
			// get profile from vall
			vall_err = false;
			core::sequence::SequenceOP vall_profile = vall.make_nmer_profile( i, frag_size, vall_err );

			// skip profiles that cross chain break or protein boundaries
			if ( !vall_err ) {
				// score query and vall profiles across all positions
				for ( Size insert_pos = 1; insert_pos <= q_prof->length() - frag_size + 1; ++insert_pos ) {
					Real score( 0.0 );
					for ( Size j = 1; j <= frag_size; ++j ) {
						Real this_score( 0.0 );
						// check to see if we know about this score
						if ( single_per_pos_scores[ insert_pos ][ i ] != -99.0 ) {
							this_score = single_per_pos_scores[ insert_pos ][ i ];
						} else {
							this_score = ss->score( q_prof, vall_profile, insert_pos + j - 1, j );
							if ( ! option[ james::debug ]() )
								single_per_pos_scores[ insert_pos ][ i ] = this_score;
						}

						score += this_score;
					}

					VallIndex idx( i, score );
					per_position_scores[ insert_pos ].add_value( idx );
				} // insert_pos
			}

			if ( i % 10000 == 0 ) {
				core::Real const percent_done( static_cast< Real > ( i / vall.size() ) );
				std::cerr << "done with " << i << " / " << vall.size() << " (" << percent_done << "%)"
									<< " scoring operations." << std::endl;
				//prof_show();
			}
		} // vall_size
		PROF_STOP( basic::SEQUENCE_COMPARISON );

		// output
		std::string output_fn =
			option[ out::file::frag_prefix ]() + "." + string_of(frag_size) + "mers";

		utility::io::ozstream output( output_fn );
		for ( core::Size insert_pos = 1;
					insert_pos <= q_prof->length() - frag_size + 1;
					++insert_pos
		) {
			output 	<<  "position: " 		<< I( 12, insert_pos )
							<< " neighbors:   " << I( 10, frags_per_position )
							<< std::endl << std::endl;
			for ( Size i = 1; i <= frags_per_position; ++i ) {
				core::Size vall_idx = per_position_scores[insert_pos][i].index();
				core::Real score    = per_position_scores[insert_pos][i].value();
				utility::vector1< std::string > lines
					= vall.make_fragment_lines( vall_idx, frag_size );

				for ( utility::vector1< std::string >::const_iterator it = lines.begin(), end = lines.end();
							it != end; ++it
				) {
					output << *it << " " << F( 8, 3, score ) << std::endl;
				} // lines
				output << std::endl;
			} // frags_per_position
		}
		output.close();
	} // frag_sizes

	std::cout << "end of program." << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
