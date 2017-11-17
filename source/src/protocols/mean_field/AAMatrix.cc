// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    AAMatrix.cc

/// @brief   Method definitions for AAMatrix.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/AAMatrix.hh>

// Package headers
#include <protocols/mean_field/RotMatrix.hh>

// Project headers
#include <core/pack/task/ResfileReader.hh>

// Utility headers
#include <protocols/mean_field/jagged_array.functions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mean_field.OptionKeys.gen.hh>

// Numeric headers


// C++ headers
#include <iostream>
#include <cmath>
#include <fstream>


namespace protocols {
namespace mean_field {


using namespace core;

static basic::Tracer TR( "protocols.mean_field.AAMatrix" );

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////

// Standard constructor
AAMatrix::AAMatrix() : jagged_array < AAProb > ()
{}

/// @details Builds AAMatrix from a standard transfac file
AAMatrix::AAMatrix( std::istream & aa_matrix_file ) : jagged_array < AAProb > ()
{
	init( aa_matrix_file );
}

/// @details Builds AAMatrix from a RotMatrix using private method init
AAMatrix::AAMatrix( RotMatrix const & rm,
	protocols::mean_field::jagged_array< core::Real > em,
	core::Real temp
): jagged_array < AAProb > ()
{
	init( rm, em, temp );
}

// Copy constructor
AAMatrix::AAMatrix( AAMatrix const & object_to_copy ) : jagged_array < AAProb > ()
{
	copy_data( *this, object_to_copy );
}

// Assignment operator
AAMatrix &
AAMatrix::operator=( AAMatrix const & object_to_copy )
{
	// Abort self-assignment.
	if ( this == &object_to_copy ) {
		return *this;
	}

	copy_data( *this, object_to_copy );
	return *this;
}

// Destructor
AAMatrix::~AAMatrix() {}

/// param [out] results - vector1 of cosine distances, calculated per-column and for the entire matrix.
/// @remarks distance is calculated between corresponding columns, determined by pos() of the AAProbs in the columns
/// @remarks if there are no corresponding columns, returns empty vector1
utility::vector1 < Real >
AAMatrix::cosine_distance( AAMatrix const & obj ) const
{
	AAMatrix pruned_curr;
	AAMatrix pruned_obj;

	for ( Size col_curr = 1; col_curr <= size(); ++col_curr ) {
		for ( Size col_obj = 1; col_obj <= obj.size(); ++col_obj ) {
			if ( ( *this )[ col_curr ][ 1 ].pos() == obj[ col_obj ][ 1 ].pos() ) {
				pruned_curr.push_back( ( *this )[ col_curr ] );
				pruned_obj.push_back( ( obj )[ col_obj ] );
			}
		}
	}

	protocols::mean_field::jagged_array < AAProb > mult_aa_matrix = pruned_curr * pruned_obj;
	protocols::mean_field::jagged_array < AAProb > curr_sq_aa_matrix = pruned_curr * pruned_curr;
	protocols::mean_field::jagged_array < AAProb > obj_sq_aa_matrix = pruned_obj * pruned_obj;

	utility::vector1 < AAProb > mult_aa_matrix_totals = mult_aa_matrix.get_totals_columns();
	utility::vector1 < AAProb > curr_sq_aa_matrix_totals = curr_sq_aa_matrix.get_totals_columns();
	utility::vector1 < AAProb > obj_sq_aa_matrix_totals = obj_sq_aa_matrix.get_totals_columns();

	mult_aa_matrix_totals.push_back( mult_aa_matrix.get_total() );
	curr_sq_aa_matrix_totals.push_back( curr_sq_aa_matrix.get_total() );
	obj_sq_aa_matrix_totals.push_back( obj_sq_aa_matrix.get_total() );

	utility::vector1 < Real > results( pruned_curr.size() + 1 );

	for ( Size ii = 1; ii <= results.size(); ++ii ) {
		Real dot_p = mult_aa_matrix_totals[ ii ].probability();
		Real magn = sqrt( curr_sq_aa_matrix_totals[ ii ].probability() ) * sqrt ( obj_sq_aa_matrix_totals[ ii ].probability() );
		results[ ii ] = dot_p / magn;
	}

	return results;

}

/// param [out] results - vector1 of Frobenius distances, calculated per-column and for the entire matrix.
/// @remarks distance is calculated between corresponding columns, determined by pos() of the AAProbs in the columns
/// @remarks if there are no corresponding columns, returns empty vector1
utility::vector1 < Real >
AAMatrix::frob_distance ( AAMatrix const & obj ) const
{
	AAMatrix pruned_curr;
	AAMatrix pruned_obj;

	for ( Size col_curr = 1; col_curr <= size(); ++col_curr ) {
		for ( Size col_obj = 1; col_obj <= obj.size(); ++col_obj ) {
			if ( ( *this )[ col_curr ][ 1 ].pos() == obj[ col_obj ][ 1 ].pos() ) {
				pruned_curr.push_back( ( *this )[ col_curr ] );
				pruned_obj.push_back( ( obj )[ col_obj ] );
			}
		}
	}

	utility::vector1 < Real > results;

	if ( ! pruned_curr.empty() ) {
		Size res_size = pruned_curr.size() + 1;

		results.assign( pruned_curr.size() + 1, Real ( 0.0 ) );

		for ( Size pos = 1; pos <= pruned_curr.size(); ++pos ) {
			for ( Size aa = 1; aa <= pruned_curr[ pos ].size(); ++aa ) {
				Real val = pow ( ( pruned_curr[ pos ][ aa ].probability() -
					pruned_obj[ pos ][ aa ].probability() ), 2.0 );
				results[ pos ] += val;
				results[ res_size ] += val;
			}
		}

		for ( Size ii = 1; ii <= results.size(); ++ii ) {
			results[ ii ] = sqrt( results[ ii ] );
		}
	}

	return results;
}

/// param [out] results - vector1 of Frobenius distances, calculated per-column and for the entire matrix.
/// @remarks distance is calculated between corresponding columns, determined by pos() of the AAProbs in the columns
/// @remarks if there are no corresponding columns, returns empty vector1
utility::vector1 < Real >
AAMatrix::ave_abs_diff ( AAMatrix const & obj ) const
{
	AAMatrix pruned_curr;
	AAMatrix pruned_obj;

	for ( Size col_curr = 1; col_curr <= size(); ++col_curr ) {
		for ( Size col_obj = 1; col_obj <= obj.size(); ++col_obj ) {
			if ( ( *this )[ col_curr ][ 1 ].pos() == obj[ col_obj ][ 1 ].pos() ) {
				pruned_curr.push_back( ( *this )[ col_curr ] );
				pruned_obj.push_back( ( obj )[ col_obj ] );
			}
		}
	}

	utility::vector1 < Real > results;

	if ( ! pruned_curr.empty() ) {
		Size res_size = pruned_curr.size() + 1;

		results.assign( res_size, Real( 0.0 ) );

		Size num_elem = num_elements( pruned_curr );

		for ( Size pos = 1; pos <= pruned_curr.size(); ++pos ) {
			for ( Size aa = 1; aa <= pruned_curr[ pos ].size(); ++aa ) {
				Real val = std::abs ( pruned_curr[ pos ][ aa ].probability() -
					pruned_obj[ pos ][ aa ].probability() );
				results[ pos ] += val / pruned_curr[ pos ].size();
				results[ res_size ] += val / num_elem;
			}
		}
	}

	return results;
}

// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods

/// @details prints representation of AAMatrix using show() function of AAProbs
/// @details in the outputted table, each column represents the probabilities of amino acids at that position
/// @details if current amino acids are set, outputs a vector of current amino acids and their probability
void
AAMatrix::show( std::ostream & output ) const
{
	output << "Specificity Profile (AAMatrix)--------------------------------------" << std::endl;

	jagged_array< AAProb >::show( output );

	utility::vector1 < AAProb > curr_aa_p = curr_aa_prob();

	bool printed_header = false;
	for ( utility::vector1< AAProb >::const_iterator iter = curr_aa_p.begin(); iter != curr_aa_p.end(); ++iter ) {
		if ( iter->pos() != Size( 0 ) ) {
			if ( ! printed_header ) {
				output << "Curr Amino Acid Probability------------------------------" << std::endl;
				printed_header = true;
			}
			output << ( *iter ) << "\t";
		}
	}
	output << std::endl;

}

/// @param [in] filename - specifies file for transfac output
/// @details output AAMatrix in transfac file format to file specified by filename
/// @details Error if filename cannot be opened
void
AAMatrix::dump_transfac( std::string const & filename ) const
{
	std::ofstream output( filename.c_str() );

	if ( ! output.is_open() ) {
		TR.Error << "Output file for dumping transfac cannot be opened " << filename << std::endl;
		return;
	}

	output << "ID Matrix" << std::endl;
	output << "PO";

	//assumes that all vectors are same size, valid assumption because they should all be size of AA:num_canonical_aas
	for ( Size aa = 1; aa <= ( *this )[ 1 ].size(); ++aa ) {
		output << "\t" << core::chemical::oneletter_code_from_aa( ( *this )[ 1 ][ aa ].aa_ind() );
	}

	for ( Size res = 1; res <= size(); ++res ) {
		output << std::endl;
		output << ( *this )[ res ][ 1 ].pos();

		for ( Size aa = 1; aa <= ( *this )[ res ].size(); ++aa ) {
			output << "\t" << ( *this )[ res ][ aa ].probability();
		}
	}
}

/// @details returns true if AAMatrix is empty of AAProbs or if AAMatrix is full of nonsense AAProbs
bool
AAMatrix::empty() const
{
	bool empty = false;
	if ( jagged_array< AAProb >::empty() ) {
		empty = true;
	} else if ( get_total().probability() == 0 ) { //AAMatrix is full of nonsense AAProbs
		empty = true;
	}
	return empty;
}


// Accessors/Mutators

/// @details returns vector with probabilities of the current amino acids based on curr_aa_ vector
utility::vector1 < AAProb >
AAMatrix::curr_aa_prob() const
{
	utility::vector1 < AAProb > ap ( size() );
	for ( Size pos = 1; pos <= curr_aa_.size(); ++pos ) {
		ap[ pos ] = ( *this )[ pos ][ curr_aa_[ pos ] ];
	}
	return ap;
}

// Private methods /////////////////////////////////////////////////////////////

// Initialize data members.
/// @details private method used to initialize AAMatrix as jagged_array of AAProbs from RotMatrix
/// @details called by constructors
/// @details also builds curr_aa_ vector
/// @param [in] RotMatrix - used to initialize AAMatrix by summing over probabilities of all rotamers for each amino acid for each position
void
AAMatrix::init( RotMatrix const & rm, protocols::mean_field::jagged_array< core::Real > em, Real temp )
{
	utility::vector1< Size > curr_rot = rm.curr_rot();

	//clears variables in order to reinitialize
	//TODO: AR: it appears that there is a method - not strictly necessary at this point, as there is no method build_aa_matrix, so that a AAMatrix cannot be reinitialized currently
	clear();
	curr_aa_.clear();

	if ( curr_rot.size() != 0 ) {
		curr_aa_.resize( rm.n_designed() );
	}

	//used to normalize post rot_norm step
	utility::vector1 < Real > totals( rm.n_designed(), Real ( 0.0 ) );

	for ( Size pos = 1; pos <= rm.size(); ++pos ) {
		//only initialize vector in AAMatrix for designable positions
		if ( rm.is_designed( pos ) ) {
			//push back vector the size of num_canonical_aas - if there is no probability of certain aa's at some positions
			//their probability will remain 0
			push_back( utility::vector1< AAProb >( core::chemical::num_canonical_aas ) );

			//position that pointing at in AAMatrix
			//may be different than pos in RotMatrix, which includes undesigned positions (only repacked)
			Size aa_matrix_pos = size();

			//if current rotamer is set
			if ( curr_aa_.size() != 0 ) {
				curr_aa_[ aa_matrix_pos ] = rm[ pos ][ curr_rot[ pos ] ].aa_ind();
			}

			//iterates through all rotamers for pos in RotMatrix, adding them to corresponding positions in AAMatrix
			for ( Size rot = 1; rot <= rm[ pos ].size(); ++rot ) {
				core::chemical::AA ind = rm[ pos ][ rot ].aa_ind();

				//hasn't been initialized to correct amino acid yet
				if ( ( *this )[ aa_matrix_pos ][ ind ].pos() == 0 ) {
					( *this )[ aa_matrix_pos ][ ind ] = AAProb ( rm[ pos ][ rot ] ); // this line is always necessary
					//     ( *this )[ aa_matrix_pos ][ ind ].probability( em[ pos ][ rot ] ); //this line is only necessary for EM
				} else {
					( *this )[ aa_matrix_pos ][ ind ] += rm[ pos ][ rot ]; // this line is only necessary for RM
					//     ( *this )[ aa_matrix_pos ][ ind ] += em[ pos ][ rot ]; // this line is only necessary for EM
					//     ( *this )[ aa_matrix_pos ][ ind ].nrot( ( *this )[ aa_matrix_pos ][ ind ].nrot() + 1); // this line is only necessary for EM

				}
			}

			//normalize for number of rotamers by dividing orig_prob/(nrot ^ rot_norm_weight)
			for ( Size aa = 1; aa <= core::chemical::num_canonical_aas; ++aa ) {
				//if AAProb wasn't set yet (i.e. AA is not allowed) set nrot, pos, and aa_ind
				if ( ( *this )[ aa_matrix_pos ][ aa ].pos() == 0 ) {
					( *this )[ aa_matrix_pos ][ aa ].pos( rm[ pos ][ 1 ].pos() );
					( *this )[ aa_matrix_pos ][ aa ].nrot( 1 ); //set to 1 to avoid division by 0
					( *this )[ aa_matrix_pos ][ aa ].aa_ind( static_cast<core::chemical::AA>( aa ) );
				}
				Size nrot = ( *this )[ aa_matrix_pos ][ aa ].nrot();
				//    nrot = nrot == 0 ? 1 : nrot; //set to 1 if this is set to 0 to avoid division by 0
				//this line is only necessary for P/nrot^RNW
				( *this )[ aa_matrix_pos ][ aa ].probability( ( *this )[ aa_matrix_pos ][ aa ].probability() / pow( nrot,
					basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] ) );
				//this line is necessary for E -> P
				//    ( *this )[ aa_matrix_pos ][ aa ].probability( exp( -1.0 * ( *this )[ aa_matrix_pos ][ aa ].probability() /
				//      temp ) );
				//this line is only necessary for nrot*RNW - E
				//    if ( ( *this )[ aa_matrix_pos ][ aa ].probability() > 100 )
				//    {
				//     ( *this )[ aa_matrix_pos ][ aa ].probability( 100 );
				//    }
				// GC expression
				//    ( *this )[ aa_matrix_pos ][ aa ].probability( exp( ( ( *this )[ aa_matrix_pos ][ aa ].nrot() *
				//      basic::options::option[ basic::options::OptionKeys::mean_field::rot_norm_weight ] -
				//      ( *this )[ aa_matrix_pos ][ aa ].probability() ) /
				//      ( temp * 100.0 ) ) );
				totals[ aa_matrix_pos ] += ( *this )[ aa_matrix_pos ][ aa ].probability();
			}
		} // if position is being designed
	} // loops through positions
	( *this ) /= totals; //rnw commented out
	TR.Debug << em[1][1] << temp << std::endl;
}

// Initialize data members.
/// @details private method used to initialize AAMatrix as jagged_array of AAProbs from input transfac file
/// @details called by constructors
/// @details curr_aa_ vector left empty
/// @param [in] aa_matrix_file - transfac file
void
AAMatrix::init( std::istream & aa_matrix_file )
{
	//clears variables in order to reinitialize
	//not strictly necessary at this point, as there is no method build_aa_matrix, so that a AAMatrix cannot be reinitialized currently
	clear();
	curr_aa_.clear();

	utility::vector1 < core::chemical::AA > aa_names( core::chemical::num_canonical_aas );
	Size lineno = 1; //first line of transfac file is ID
	utility::vector1 < std::string > tokens( core::pack::task::tokenize_line ( aa_matrix_file) );

	while ( !tokens.empty() ) {

		if ( lineno == 2 ) { //second line of transfac file holds AA codes
			parse_aa_line( tokens, aa_names );
		} else if ( lineno > 2 ) { //subsequent lines consist of a position and probabilities for all amino acids for that position
			push_back( parse_aa_matrix_line( tokens, aa_names ) );
		}
		tokens = core::pack::task::tokenize_line ( aa_matrix_file) ;
		++lineno;

	}
}

/// @details parses line of transfac i.e. 200 0.5 0.2 0.3...
/// @details pushes back appropriate vector into AAMatrix
//TODO: AR: more error checking (limits 0<=prob<=1, length of each line is the same, etc.)
utility::vector1 < AAProb >
AAMatrix::parse_aa_matrix_line( utility::vector1 < std::string > const & tokens,
	utility::vector1 < core::chemical::AA > const & aa_names)
{

	std::istringstream token_s;
	Size pos;
	token_s.str(tokens[1]);

	if ( ! ( token_s >> pos ) ) {
		std::stringstream error_message;
		error_message
			<< "Error parsing specificity profile (AAMatrix): expected numeric value." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
	}

	utility::vector1 < AAProb > probs( core::chemical::num_canonical_aas );

	//tokens begins with positional index
	for ( Size which_token = 2; which_token <= tokens.size(); ++which_token ) {
		std::istringstream token_s( tokens[ which_token ] );
		Real prob;

		if ( ! ( token_s >> prob ) ) {
			std::stringstream error_message;
			error_message
				<< "Error parsing specificity profile (AAMatrix): expected numeric value." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message.str() );
		}

		core::chemical::AA aa = aa_names[ which_token - 1 ];

		probs[aa] = AAProb( prob, aa, pos, 0 );
	}

	return probs;
}

/// @details parses header (2nd line) of transfac file into list of aa_names
/// @details necessary to determine which probability is which AA for subsequent lines
void
AAMatrix::parse_aa_line( utility::vector1< std::string > const & tokens, utility::vector1< core::chemical::AA > & aa_names )
{
	for ( Size which_token = 2; which_token <= tokens.size(); ++which_token ) {
		aa_names[ which_token - 1 ] = core::chemical::aa_from_oneletter_code( tokens[ which_token ].at( 0 ) );
	}
}


// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
AAMatrix::copy_data(
	AAMatrix & object_to_copy_to,
	AAMatrix const & object_to_copy_from)
{
	object_to_copy_to.clear();
	for ( Size i = 1; i <= object_to_copy_from.size(); ++i ) {
		object_to_copy_to.push_back ( object_to_copy_from[i] );
	}
	object_to_copy_to.curr_aa_ = object_to_copy_from.curr_aa_;
}


// Friend methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that AAMatrix can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, AAMatrix const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

} //namespace mean_field
}  // namespace protocols
