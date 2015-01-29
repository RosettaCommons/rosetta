// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SequenceProfile.cc
/// @brief class definition for a given scoring scheme for an alignment.
/// @detailed Simply based on comparing single characters from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/MatrixScoringScheme.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/chemical/AA.hh>

// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <iostream>
#include <string>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

namespace core {
namespace sequence {

static thread_local basic::Tracer tr( "core.sequence.SequenceProfile" );

void SequenceProfile::read_from_checkpoint(
		utility::file::FileName const & fn,
		bool negative_better
) {
	utility::io::izstream input( fn );
	negative_better_ = negative_better;
	// order of amino acids in the .checkpoint file
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

	std::string aa_seq;
	utility::vector1< utility::vector1< core::Real > > new_prof;
	// profile is indexed by the order of amino acids in core::chemical::AA

	if ( !input ) {
		utility_exit_with_message( "ERROR: Unable to open file!" );
	}
	std::string line;

	//Check format of file header
	getline( input, line );

	utility::vector1< std::string > header( utility::split_whitespace( line ) );
	if( header.size() == order.size() || header.size() == order.size()+1 ) { // Skip header if it doesn't match the expected number of columns.
		core::Size offset(0);
		if( header.size() == order.size()+1 ) { offset = 1; }
		for( core::Size ii=1+offset; ii <= header.size(); ++ii) {
			if( header[ii].size() != 1 ) { break; } // Not a single letter code
			if( core::chemical::oneletter_code_from_aa(order[ii-offset]) != header[ii][0] ) {
				tr.Warning << "WARNING: Potential badly formatted sequence profile file '" << std::string(fn) << "'. " <<
						"Columns should be in ACDEFGHIKLMNPQRSTVWY order. " <<
						"Saw '" << header[ii] << "' as column header at position " << utility::to_string(ii) << ". " <<
						"Expected '" << core::chemical::oneletter_code_from_aa(order[ii-offset]) << "'." << std::endl;
				break;
			}
			alphabet_.push_back( header[ii] );
		}
	}
	if( alphabet_.size() != order.size() ) {
			tr.Debug << "No header for sequence profile checkpoint file '" << fn << "'. Assuming ACDEFGHIKLMNPQRSTVWY order." << std::endl;
	}

	//Get data in body of file
	while( getline( input, line ) ) {
		if ( line.substr(0,3) == "END" ) break; // end of profile
		std::string aa;
		utility::vector1< core::Real > prof_row;
		prof_row.resize( order.size() );

		std::istringstream ls( line );
		ls >> aa;

		core::Real aa_prob;
		ls >> aa_prob;
		Size index(1);
		while ( !ls.fail() ) {
			//prof_row.push_back( aa_prob );
			prof_row[ order[index] ] = aa_prob;

			ls >> aa_prob;
			++index;
		}

		if( index != order.size() + 1 ) {
			tr.Warning << "WARNING: Potentially incomplete profile row in '" << std::string(fn) << "'. " << std::endl;
		}
		aa_seq += aa;
		new_prof.push_back( prof_row );
	}

	sequence( aa_seq );
	profile( new_prof );
}

void SequenceProfile::read_from_binary_chk(utility::file::FileName const & fn) {

    double x = 0;
    char bb4[4];
    std::ifstream myfile (fn.name().c_str(), std::ios::in|std::ios::binary);
    if ( !myfile ) {
		std::string msg(
			"ERROR: Unable to open file " +
			static_cast< std::string > (fn) +
			"!"
		);
		utility_exit_with_message( msg );
    }

	static utility::vector1< core::chemical::AA > order;
	order.resize( 20 );
	order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
	order[ 2] = core::chemical::aa_from_oneletter_code( 'R' );
	order[ 3] = core::chemical::aa_from_oneletter_code( 'N' );
	order[ 4] = core::chemical::aa_from_oneletter_code( 'D' );
	order[ 5] = core::chemical::aa_from_oneletter_code( 'C' );
	order[ 6] = core::chemical::aa_from_oneletter_code( 'Q' );
	order[ 7] = core::chemical::aa_from_oneletter_code( 'E' );
	order[ 8] = core::chemical::aa_from_oneletter_code( 'G' );
	order[ 9] = core::chemical::aa_from_oneletter_code( 'H' );
	order[10] = core::chemical::aa_from_oneletter_code( 'I' );
	order[11] = core::chemical::aa_from_oneletter_code( 'L' );
	order[12] = core::chemical::aa_from_oneletter_code( 'K' );
	order[13] = core::chemical::aa_from_oneletter_code( 'M' );
	order[14] = core::chemical::aa_from_oneletter_code( 'F' );
	order[15] = core::chemical::aa_from_oneletter_code( 'P' );
	order[16] = core::chemical::aa_from_oneletter_code( 'S' );
	order[17] = core::chemical::aa_from_oneletter_code( 'T' );
	order[18] = core::chemical::aa_from_oneletter_code( 'W' );
	order[19] = core::chemical::aa_from_oneletter_code( 'Y' );
	order[20] = core::chemical::aa_from_oneletter_code( 'V' );


    myfile.read(bb4,4);

    int seqLength = bb4[0];
    int b2 = ((int) bb4[1]) << 8;
    int b3 = ((int) bb4[2]) << 16;
    int b4 = ((int) bb4[3]) << 24;
    seqLength = (seqLength < 0) ? 256 + seqLength + b2 + b3 + b4 : seqLength + b2 + b3 + b4;
    std::cout << seqLength<<" "<<(int)bb4[0]<<" "<<(int)bb4[1] << "\n";

    utility::vector1<char> seq;
    for(int i=0;i<seqLength;++i) {
      seq.push_back( myfile.get() );
    }
    string strSeq(seq.begin(),seq.end());
    tr.Debug << "Read sequence " << strSeq << " from  " << fn << std::endl;

    char bb8[8];
    //double row[20]; /// @ralford werror catches as an unused var 3/26/14
    for(int i = 0;i < seqLength;i++) {
      utility::vector1< Real > prof_row;
      for(int j = 0;j < 20;j++) {
        myfile.read(bb8,8);
        std::copy(bb8, bb8 + sizeof(double), reinterpret_cast<char*>(&x));
//        x = swap(x);
        std::cout << std::fixed<<std::setw(7)<<std::setprecision(4)<<x<<" ";
        prof_row.push_back(x);
      }
      profile_.push_back( prof_row );
      std::cout<<std::endl;
    }
}

void SequenceProfile::read_from_file(
	utility::file::FileName const & fn
) {
	negative_better_ = false;
	// order of amino acids in the .pssm file
	static utility::vector1< core::chemical::AA > order;
	order.resize( 20 );
	order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
	order[ 2] = core::chemical::aa_from_oneletter_code( 'R' );
	order[ 3] = core::chemical::aa_from_oneletter_code( 'N' );
	order[ 4] = core::chemical::aa_from_oneletter_code( 'D' );
	order[ 5] = core::chemical::aa_from_oneletter_code( 'C' );
	order[ 6] = core::chemical::aa_from_oneletter_code( 'Q' );
	order[ 7] = core::chemical::aa_from_oneletter_code( 'E' );
	order[ 8] = core::chemical::aa_from_oneletter_code( 'G' );
	order[ 9] = core::chemical::aa_from_oneletter_code( 'H' );
	order[10] = core::chemical::aa_from_oneletter_code( 'I' );
	order[11] = core::chemical::aa_from_oneletter_code( 'L' );
	order[12] = core::chemical::aa_from_oneletter_code( 'K' );
	order[13] = core::chemical::aa_from_oneletter_code( 'M' );
	order[14] = core::chemical::aa_from_oneletter_code( 'F' );
	order[15] = core::chemical::aa_from_oneletter_code( 'P' );
	order[16] = core::chemical::aa_from_oneletter_code( 'S' );
	order[17] = core::chemical::aa_from_oneletter_code( 'T' );
	order[18] = core::chemical::aa_from_oneletter_code( 'W' );
	order[19] = core::chemical::aa_from_oneletter_code( 'Y' );
	order[20] = core::chemical::aa_from_oneletter_code( 'V' );

	utility::io::izstream input( fn );

	if ( !input ) {
		std::string msg(
			"ERROR: Unable to open file " +
			static_cast< std::string > (fn) +
			"!"
		);
		utility_exit_with_message( msg );
	}
	std::string line;

	tr.Debug << "reading from " << fn << std::endl;
	//tr << "reading from " << fn << std::endl;
	// read in two header lines
	getline( input, line );
	getline( input, line );

	// initialize headers
	getline( input, line );
	std::istringstream line_stream( line );
	while ( !line_stream.fail() ) {
		std::string aa;
		line_stream >> aa;
		if ( line_stream.fail() ) continue;

		alphabet_.push_back( aa );

		if( core::chemical::oneletter_code_from_aa(order[alphabet_.size()]) != aa[0] ) {
			utility_exit_with_message("Badly formatted pssm file. Expected ARNDCQEGHILKMFPSTWYV order in third line of file '"+
					std::string(fn)+"'. Saw '"+aa+"' at position "+utility::to_string(alphabet_.size())+". "
					"Expected '"+core::chemical::oneletter_code_from_aa(order[alphabet_.size()])+"'.");
		}
		if ( alphabet_.size() >= 20 ) break; // super-hack for pssm file format! bad james, bad!
	}

	std::string seq;
	while( getline( input, line ) ) {
		std::istringstream line_stream( line );
		//std::cout << "line = " << line << std::endl;
		Size pos;
		string aa;

		line_stream >> pos >> aa;
	//	tr <<"pos: "<<pos<<"aa: "<<aa<<std::endl;
		if ( line_stream.fail() ) continue;

		utility::vector1< Real > prof_row;
		utility::vector1< Real > probability_row; //will hold the probabilty vector, gideonla 120214
		prof_row.resize( order.size() );
		probability_row.resize( order.size() );


		Real score;
		line_stream >> score;
		Size index(1);
		while ( !line_stream.fail() && index <= order.size() ) {
			prof_row[ order[index] ] = score;
			line_stream >> score;
			++index;
		}


		//using the same while loop to read % matrix in the pssm file. All the numbers after index 20 are stored in the percentage matrix. I am doing it this way
		// So I don't have to change anything in the existing code , gideonla 120214
		if (basic::options::option[ basic::options::OptionKeys::out::file::use_occurrence_data ].value()){
		index=1;
		while ( !line_stream.fail() && index <= order.size() ) {
		//	tr <<"Percentage :"<<score<<std::endl;
			probability_row[ order[index] ] = score;
			line_stream >> score;
			++index;
		}// end of the while loop
		if (index<order.size()){//if the probability matrix is missing or corrupt we fill the the % matrix with zeros
			for (core::Size i(1), end(order.size()); i <= end; ++i) {
				probability_row[order[i]] = 0;
			}
		}
		occurrence_data_.push_back( probability_row );
		}

		profile_.push_back( prof_row );


		seq += aa;
	} // while( getline( input, line ) )

	if( profile_.size() == 0 ) {
		utility_exit_with_message("Profile file '"+std::string(fn)+"' does not appear to contain any data.");
	}

	tr.Debug << "Read sequence " << seq << " from  " << fn << std::endl;
	tr.Debug << "profile dimensions are " << profile_.size() << "x"
		<< profile_.front().size() << "." << std::endl;
	sequence( seq );
	id( std::string(fn) );
	//tr<<" the size of profile_ is: "<<profile_.size()<<" and the size of occurrence_data_ is: "<<occurrence_data_.size()<<std::endl;
} // read_from_file



void SequenceProfile::generate_from_sequence( Sequence const & seq, std::string matrix ) {
	MatrixScoringScheme mat;
	mat.read_from_database(matrix);

	utility::vector1< utility::vector1< core::Real > > new_prof;
	for (core::Size ii(1), end(seq.length()); ii <= end; ++ii) {
		utility::vector1< Real > prof_row( mat.values_for_aa( seq[ii] ) );
		prof_row.resize( core::chemical::num_canonical_aas, -10000 );
		new_prof.push_back(prof_row);
	}

	sequence( seq.sequence() );
	profile( new_prof );
	negative_better_ = false;
}

/// @brief Returns the 2D vector1 of Real-values representing this profile.
utility::vector1< utility::vector1< Real > > const &
SequenceProfile::profile() const {
	return profile_;
}

utility::vector1< utility::vector1< Real > > const &
SequenceProfile::occurrence_data() const {
	return occurrence_data_;
}

/// @brief Multiply all profile weights by factor
void SequenceProfile::rescale(core::Real factor) {
	for ( Size ii = 1; ii <= profile().size(); ++ii ) {
		for ( utility::vector1< core::Real >::iterator
						it = profile_[ii].begin(), end = profile_[ii].end(); it != end; ++it ) {
			*it *= factor;
		}
	}
	if ( factor < 0 ) {
		tr << "Flipping sense of negative_better." << std::endl;
		negative_better_ = ! negative_better_;
	}
}

void SequenceProfile::convert_profile_to_probs( core::Real temp /*= 1.0*/ ) {
	temp_ = temp;
	utility::vector1< utility::vector1< Real > > new_prof;
	for ( Size ii = 1; ii <= profile().size(); ++ii ) {
		utility::vector1< core::Real > new_prof_row = prof_row(ii);
		scores_to_probs_( new_prof_row, temp, negative_better_ );
		new_prof.push_back( new_prof_row );
	}
	profile( new_prof );
	negative_better_ = false;
}

void SequenceProfile::global_auto_rescale() {
	core::Real maxval(0.0);
	for ( Size ii = 1; ii <= profile_.size(); ++ii ) {
		for ( Size jj = 1; jj <= profile_[ii].size(); ++jj ) {
			if ( maxval < std::abs(profile_[ii][jj]) ) {
				maxval = std::abs(profile_[ii][jj]);
			}
		}
	}
	rescale( 1/maxval );
}


void SequenceProfile::profile(
	utility::vector1< utility::vector1< core::Real > > const & new_prof
) {
	profile_ = new_prof;
}
void SequenceProfile::occurrence_data(
	utility::vector1< utility::vector1< core::Real > > const & new_occurrence_data
) {
	occurrence_data_ = new_occurrence_data;
}

void SequenceProfile::prof_row( utility::vector1< Real > const & new_prof_row, core::Size pos ) {
	if( profile_.size() <= pos )
		profile_.resize(pos+1);
	profile_[pos] = new_prof_row;
}

void SequenceProfile::probabilty_row( utility::vector1< Real > const & new_prob_row, core::Size pos ) {
	if( occurrence_data_.size() <= pos )
		occurrence_data_.resize(pos+1);
	occurrence_data_[pos] = new_prob_row;
}

void SequenceProfile::insert_char(
	core::Size pos,
	char new_char
) {
	using core::Real;
	using utility::vector1;
	// add in a profile column of zeroes.
	vector1< Real > zero_col( width(), 0.0 );
	vector1< vector1< Real > > new_prof, old_prof( profile() );
	for ( Size i = 0; i <= length() + 1; ++i ) {
		if ( i == pos )                new_prof.push_back( zero_col );
		if ( i >= 1 && i <= length() ) new_prof.push_back( old_prof[i] );
	}

	profile( new_prof );
	Sequence::insert_char( pos, new_char );
}

void SequenceProfile::delete_position(
	core::Size pos
) {
	//std::cout << "sequence is " << sequence() << std::endl;
	using core::Real;
	using utility::vector1;

	runtime_assert( pos <= length() );

	vector1< vector1< Real > > new_prof( profile() );
	vector1< vector1< Real > >::iterator it
		= new_prof.begin() + pos - start();

	runtime_assert( it != new_prof.end() );

	new_prof.erase( it );
	profile( new_prof );

	Sequence::delete_position( pos );
	//std::cout << "sequence is " << sequence() << std::endl;
}

/// @brief Returns the number of distinct values at each position in this profile.
Size SequenceProfile::width() const {
debug_assert( check_internals_() );
	return profile_[1].size();
}

/// @brief Returns the vector1 of values at this position.
utility::vector1< Real > const &
SequenceProfile::prof_row( Size pos ) const {
	runtime_assert( pos <= profile_.size() );
	return profile_[ pos ];
}

/// @brief Returns the vector1 of values at this position.
utility::vector1< Real > const &
SequenceProfile::probability_row( Size pos ) const {
	runtime_assert( pos <= occurrence_data_.size() );
	return occurrence_data_[ pos ];
}

void SequenceProfile::scores_to_probs_(
	utility::vector1< core::Real > & scores,
	core::Real kT,
	bool negative_better /* = false */
) const {
	using std::exp;
	using utility::vector1;

	// calculate partition (aka Z), with this definition:
	// Z = sum( exp( score / kT ) )
	core::Real partition( 0.0 );
	for ( vector1< core::Real >::iterator
				it = scores.begin(), end = scores.end(); it != end; ++it
	) {
		if ( negative_better ) {
			*it = exp( -1 * *it / kT );
		} else {
			*it = exp( *it / kT );
		}
		partition += *it;
	}

	// transform scores using the partition calculated above:
	// P(s) = exp( -1  * score / kT ) ) / Z
	for ( utility::vector1< core::Real >::iterator
				it = scores.begin(), end = scores.end(); it != end; ++it
	) {
		*it = *it / partition;
	}
} // scores_to_probs_

std::ostream & operator<<( std::ostream & out, const SequenceProfile & p ) {
	Size width = 8;
	Size precision = 3;

	out << p.to_string() << std::endl;
	for ( Size i = 1; i <= p.length(); ++i ) {
		for ( Size j = 1; j <= p.width(); ++j ) {
			out << ObjexxFCL::format::F( width, precision, p.prof_row(i)[j] );
		}
		out << std::endl;
	}

	return out;
}

bool SequenceProfile::check_internals_() const {
	using core::Real;
	using utility::vector1;

	runtime_assert( profile_.size() == length() );

	for ( Size i = 1; i <= profile_.size(); ++i ) {
		runtime_assert( profile_[i].size() == alphabet_.size() );
	}

	return true;
}

} // sequence
} // core
