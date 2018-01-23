// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.cc
/// @brief A helper class to assist in parsing rotamer libraries.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

//Core includes
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.tmpl.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/chemical/ResidueType.hh>

//Basic includes
#include <basic/Tracer.hh>

//C++ includes
#include <iostream>

static basic::Tracer TR( "core.pack.dunbrack.RotamericSingleResidueDunbrackLibraryParser" );


namespace core {
namespace pack {
namespace dunbrack {

#define MIN_CHI_STD_DEV 0 //Below this, the chi standard deviation is replaced with the bogus value.
#define BOGUS_CHI_STD_DEV 5 //The bogus value used if chi_std_dev is too small.
#define MIN_PROBABILITY 1e-6 //The minimum rotamer probability.  Below this value, the value is set to this value.

/// @brief Default constructor.
/// @param[in] num_mainchain_torsions The number of mainchain torsions on which this rotamer library depends.
/// @param[in] max_possible_chis The maximum number of chis allowed in this rotamer library.  This is the number of chi columns in the file (usually 4).
/// @param[in] max_possible_rotamers The maximum number of rotamers that this rotamer library can define.
RotamericSingleResidueDunbrackLibraryParser::RotamericSingleResidueDunbrackLibraryParser(
	core::Size const num_mainchain_torsions,
	core::Size const max_possible_chis,
	core::Size const max_possible_rotamers
):
	utility::pointer::ReferenceCount(),
	read_file_was_called_(false),
	num_mainchain_torsions_(num_mainchain_torsions),
	max_possible_chis_(max_possible_chis),
	max_possible_rotamers_(max_possible_rotamers),
	backbone_torsions_(),
	counts_(),
	rotwells_(),
	probabilities_(),
	chimeans_(),
	chi_std_devs_()
{
	debug_assert( num_mainchain_torsions_ > 0 );
	debug_assert( max_possible_chis_ > 0);
	debug_assert( max_possible_rotamers_ > 0 );
}

RotamericSingleResidueDunbrackLibraryParser::~RotamericSingleResidueDunbrackLibraryParser(){}

RotamericSingleResidueDunbrackLibraryParser::RotamericSingleResidueDunbrackLibraryParser( RotamericSingleResidueDunbrackLibraryParser const &src ) :
	utility::pointer::ReferenceCount(),
	read_file_was_called_(src.read_file_was_called_),
	num_mainchain_torsions_(src.num_mainchain_torsions_),
	max_possible_chis_(src.max_possible_chis_),
	max_possible_rotamers_(src.max_possible_rotamers_),
	backbone_torsions_(src.backbone_torsions_),
	counts_(src.counts_),
	rotwells_(src.rotwells_),
	probabilities_(src.probabilities_),
	chimeans_(src.chimeans_),
	chi_std_devs_(src.chi_std_devs_)
{
	debug_assert( num_mainchain_torsions_ > 0 );
	debug_assert( max_possible_chis_ > 0);
	debug_assert( max_possible_rotamers_ > 0);
}


RotamericSingleResidueDunbrackLibraryParserOP
RotamericSingleResidueDunbrackLibraryParser::clone() const {
	return RotamericSingleResidueDunbrackLibraryParserOP( new RotamericSingleResidueDunbrackLibraryParser( *this ) );
}

/// @brief Clear all of the data read from a file.
/// @details Also resets read_file_was_called_ to false.
void
RotamericSingleResidueDunbrackLibraryParser::clear_imporated_data() {
	read_file_was_called_ = false;

	backbone_torsions_.clear();
	counts_.clear();
	rotwells_.clear();
	probabilities_.clear();
	chimeans_.clear();
	chi_std_devs_.clear();

	backbone_torsions_.reserve(max_possible_rotamers_);
	counts_.reserve(max_possible_rotamers_);
	rotwells_.reserve(max_possible_rotamers_);
	probabilities_.reserve(max_possible_rotamers_);
	chimeans_.reserve(max_possible_rotamers_);
	chi_std_devs_.reserve(max_possible_rotamers_);
}

/// @brief Read an input stream to set up this object.
/// @details The RotamericSingleResidueDunbrackLibraryParser stores the information from the rotamer library in a very one-to-one
/// format, and can subsequently be used to configure a RotamericSingleResidueDunbrackLibrary object (requiring some interpretation
/// of the data).
/// @note We might be starting the read in the middle of the file, and the three-letter code from the first line might already have been read.
/// @returns The name of the next amino acid specified in the input stream, or alternatively an empty string if EOF is reached.
std::string
RotamericSingleResidueDunbrackLibraryParser::read_file(
	utility::io::izstream & infile,
	bool first_line_three_letter_code_already_read,
	core::chemical::AA const aa
) {
	//Clear all existing data:
	clear_imporated_data();

	//Indicate that we have now called this function:
	read_file_was_called_ = true;

	bool first_noncomment_line(true); //Switched to false after the first line that isn't a comment is read.

	std::string const my_name( chemical::name_from_aa( aa ) );
	std::string three_letter_code;
	bool is_shapovalov_file(false);
	core::Size read_rotamer_count(0); //Number of rotamers successfully read in.

	while ( infile ) {

		/// 1a. peek at the line; if it starts with #, skip to the next line.
		/// Also, set up an istringstream for the line.
		char first_char = infile.peek();
		if ( first_char == '#' ) {
			std::string line;
			infile.getline( line );
			first_line_three_letter_code_already_read = false;
			continue;
		}

		/// 1b. Shapovalov files have an extra column in the MIDDLE, which is a pain in the neck.
		/// First, we need to determine whether that applies here.
		/// (Added by VKM, 11 July 2016).
		if ( first_noncomment_line ) {
			first_noncomment_line = false; //We only do this once per file, since it requires two parses of the same line.  We use the first line that isn't commented out.
			is_shapovalov_file = check_for_extra_column(infile, first_line_three_letter_code_already_read);
		}

		// 2. Read the line.  Format is:
		// a.   three-letter-code
		// b.   phi
		// c.   psi (note that there can be more than two mainchain torsions, here)
		// d.   count
		// e..h r1..r4 (only used for warnings)
		// i.   probability
		// iprime. neg ln probability (iff Shapovalov)
		// j..m chimean1..chimean4
		// n..q chisd1..chisd4

		// a.  Read the three-letter code, if we haven't already.
		if ( first_line_three_letter_code_already_read ) {
			/// The last library to read from this file already ate my three letter code;
			/// skip directly to reading phi/psi values...
			first_line_three_letter_code_already_read = false;
		} else {
			infile >> three_letter_code;
			// AMW: this is an issue. What if we are a variant lib like TRP,
			// but our TLC is UNK? thankfully we can add in that exact condition.
			if ( infile.eof() || ( three_letter_code != my_name && three_letter_code != "UNK" ) ) { //If this is true, then we're done reading the current amino acid, so return the name of the next.
				check_correct_vector_lengths();
				if ( infile.eof() ) return "";
				return three_letter_code;
			} // else, we're still reading data for the intended amino acid, so let's continue...
		}

		// b, c.  Read the backbone torsions.
		{
			utility::vector1< core::Real > bb( num_mainchain_torsions_ );
			for ( Size ii = 1; ii <= num_mainchain_torsions_; ++ii ) {
				infile >> bb[ ii ];
			}
			backbone_torsions_.push_back(bb);
		}

		// d,e-h,i.  Read count, rotwells, and probability.
		{
			core::Size count;
			infile >> count;
			counts_.push_back(count);

			utility::vector1< core::Size > rotwells( max_possible_chis_ );
			for ( core::Size i(1); i<=max_possible_chis_; ++i ) {
				infile >> rotwells[ i ];
			}
			rotwells_.push_back(rotwells);

			core::Real probability;
			infile >> probability;
			if ( probability <= MIN_PROBABILITY ) probability = MIN_PROBABILITY;
			// APL -- On the advice of Roland Dunbrack, modifying the minimum probability to the
			// resolution of the library.  This helps avoid overwhelmingly unfavorable energies
			// (5 log-units difference between 1e-4 and 1e-9) for rare rotamers.
			// AMW -- Changing this to be a <= criterion and to 1e-6 to actually match the resolution
			// of the library
			probabilities_.push_back(probability);
		}

		// iprime.  Read -ln(prob), if this is a Shapovalov file.
		if ( is_shapovalov_file ) {
			core::Real minuslnProb;
			infile >> minuslnProb; //Not actually used for anything, but must be parsed.  Irritating.
		}

		// j-m.  The chimeans.  These are actually important.
		{
			utility::vector1< core::Real > chimean( max_possible_chis_ );
			for ( core::Size i(1); i<=max_possible_chis_; ++i ) {
				infile >> chimean[i];
			}
			chimeans_.push_back(chimean);
		}

		// n-q.  The chistandard deviations.  These are also important.
		{
			utility::vector1< core::Real > chi_std_devs( max_possible_chis_ );
			for ( core::Size i(1); i<=max_possible_chis_; ++i ) {
				infile >> chi_std_devs[i];
				if ( chi_std_devs[i] <= MIN_CHI_STD_DEV ) chi_std_devs[i] = BOGUS_CHI_STD_DEV; //Avoid inf.
			}
			chi_std_devs_.push_back(chi_std_devs);
		}

		if ( !correct_for_duplicated_mainchain_data() ) ++read_rotamer_count;

	} //Continue to next line in file.

	check_correct_vector_lengths();

	return ""; //If we've arrived here, then we reached an EOF, so return an empty string for the name of the next amino acid in the input stream (i.e. there is no next one).
}

/////PRIVATE MEMBER FUNCTIONS/////

/// @brief When the first non-comment line is read from a rotamer file, check whether there's an extra column (signifying that this is a
/// Shapovalov file, which has an extra column in the middle.)
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
bool
RotamericSingleResidueDunbrackLibraryParser::check_for_extra_column(
	utility::io::izstream & infile,
	bool const first_line_three_letter_code_already_read
) const {
	std::string line;
	infile.getline( line );
	//This is SUPER hacky, but need to rewind to the start of the line, and without using seekg or tellg (which aren't available for gzipped streams):
	for ( core::Size linelength( line.length() + 1 ) ; linelength > 0; --linelength ) {
		infile.stream().unget();
	}

	//Need to remove whitespace:
	utility::trim( line, " \t" );
	core::Size const shift( first_line_three_letter_code_already_read ? 1 : 0 );

	std::istringstream curline(line);
	std::string dummy; //Dummy string as recepticle for parsing line.
	core::Size counter(0);
	//Count how many whitespace-separated things there are in the line:
	while ( !curline.eof() ) {
		curline >> dummy;
		++counter;
	}
	if ( counter == 11 + num_mainchain_torsions_ + max_possible_chis_ - shift ) return false;
	runtime_assert_string_msg( counter == 12 + num_mainchain_torsions_ + max_possible_chis_ - shift, "Error in core::pack::dunbrack::RotamericSingleResidueDunbrackLibraryParser::check_for_extra_column(): A rotamer file must have either 15+N or 16+N columns (where N is the number of mainchain torsions upon which the rotamer library depends)." );
	return true;
}

/// @brief The legacy reader had a check for mainchain torsions equal to -180, in which case the last rotamer read was discarded.
/// Maintaining that here.
/// @returns Returns true for correction applied (i.e. last rotamer deleted), false otherwise.
bool
RotamericSingleResidueDunbrackLibraryParser::correct_for_duplicated_mainchain_data() {
	utility::vector1< core::Real > const & last_bb( backbone_torsions_[backbone_torsions_.size()] );
	debug_assert(last_bb.size() == num_mainchain_torsions_); //Should be true
	bool duplicated(false);
	for ( core::Size i(1); i<=num_mainchain_torsions_; ++i ) {
		if ( last_bb[i] == -180 /*exact equality -- mildly risky*/ ) {
			duplicated=true;
			break;
		}
	}
	if ( duplicated ) {
		remove_last_entry_in_vectors();
	}
	return duplicated;
}

/// @brief Delete the last entry in all of the data vectors.
void
RotamericSingleResidueDunbrackLibraryParser::remove_last_entry_in_vectors() {
	backbone_torsions_.pop_back();
	counts_.pop_back();
	rotwells_.pop_back();
	probabilities_.pop_back();
	chimeans_.pop_back();
	chi_std_devs_.pop_back();
}

/// @brief Confirm that all data vectors have the same length.  Throw an error if they don't.
void
RotamericSingleResidueDunbrackLibraryParser::check_correct_vector_lengths() const {
	core::Size const rightsize( backbone_torsions_.size() );
	static const std::string errmsg( "Error in core::pack::dunbrack::RotamericSingleResidueDunbrackLibraryParser::check_correct_vector_lengths(): The " );
	static const std::string errmsg2( " vector is the wrong length." );
	runtime_assert_string_msg( counts_.size() == rightsize, errmsg + "counts_" + errmsg2 );
	runtime_assert_string_msg( rotwells_.size() == rightsize, errmsg + "rotwells_ " + errmsg2 );
	runtime_assert_string_msg( probabilities_.size() == rightsize, errmsg + "probabilities_ " + errmsg2 );
	runtime_assert_string_msg( chimeans_.size() == rightsize, errmsg + "chimeans_ " + errmsg2 );
	runtime_assert_string_msg( chi_std_devs_.size() == rightsize, errmsg + "chi_std_devs_ " + errmsg2 );
}

/// @brief THIS FUNCTION SHOULD NEVER BE CALLED!
void
RotamericSingleResidueDunbrackLibraryParser::hacky_parser_init_workaround() {
	utility_exit_with_message( "The hacky_parser_init_workaround() function should never be called!" );

	core::chemical::ResidueType rt( nullptr, nullptr, nullptr, nullptr );

	RotamericSingleResidueDunbrackLibrary< 1, 1 > lib_1_1( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 1, 2 > lib_1_2( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 1, 3 > lib_1_3( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 1, 4 > lib_1_4( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 1, 5 > lib_1_5( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 2, 1 > lib_2_1( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 2, 2 > lib_2_2( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 2, 3 > lib_2_3( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 2, 4 > lib_2_4( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 2, 5 > lib_2_5( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 3, 1 > lib_3_1( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 3, 2 > lib_3_2( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 3, 3 > lib_3_3( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 3, 4 > lib_3_4( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 3, 5 > lib_3_5( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 4, 1 > lib_4_1( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 4, 2 > lib_4_2( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 4, 3 > lib_4_3( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 4, 4 > lib_4_4( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 4, 5 > lib_4_5( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 5, 1 > lib_5_1( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 5, 2 > lib_5_2( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 5, 3 > lib_5_3( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 5, 4 > lib_5_4( rt, false, false, false, 0.1, 0.1, false );
	RotamericSingleResidueDunbrackLibrary< 5, 5 > lib_5_5( rt, false, false, false, 0.1, 0.1, false );

	RotamericSingleResidueDunbrackLibraryParser parser_1_1( 1, 1, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_1_2( 1, 2, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_1_3( 1, 3, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_1_4( 1, 4, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_1_5( 1, 5, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_2_1( 2, 1, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_2_2( 2, 2, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_2_3( 2, 3, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_2_4( 2, 4, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_2_5( 2, 5, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_3_1( 3, 1, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_3_2( 3, 2, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_3_3( 3, 3, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_3_4( 3, 4, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_3_5( 3, 5, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_4_1( 4, 1, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_4_2( 4, 2, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_4_3( 4, 3, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_4_4( 4, 4, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_4_5( 4, 5, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_5_1( 5, 1, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_5_2( 5, 2, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_5_3( 5, 3, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_5_4( 5, 4, 1 );
	RotamericSingleResidueDunbrackLibraryParser parser_5_5( 5, 5, 1 );

	parser_1_1.configure_rotameric_single_residue_dunbrack_library< 1, 1 >( lib_1_1, utility::fixedsizearray1< core::Size, 1 >(0) );
	parser_1_2.configure_rotameric_single_residue_dunbrack_library< 1, 2 >( lib_1_2, utility::fixedsizearray1< core::Size, 2 >(0) );
	parser_1_3.configure_rotameric_single_residue_dunbrack_library< 1, 3 >( lib_1_3, utility::fixedsizearray1< core::Size, 3 >(0) );
	parser_1_4.configure_rotameric_single_residue_dunbrack_library< 1, 4 >( lib_1_4, utility::fixedsizearray1< core::Size, 4 >(0) );
	parser_1_5.configure_rotameric_single_residue_dunbrack_library< 1, 5 >( lib_1_5, utility::fixedsizearray1< core::Size, 5 >(0) );
	parser_2_1.configure_rotameric_single_residue_dunbrack_library< 2, 1 >( lib_2_1, utility::fixedsizearray1< core::Size, 1 >(0) );
	parser_2_2.configure_rotameric_single_residue_dunbrack_library< 2, 2 >( lib_2_2, utility::fixedsizearray1< core::Size, 2 >(0) );
	parser_2_3.configure_rotameric_single_residue_dunbrack_library< 2, 3 >( lib_2_3, utility::fixedsizearray1< core::Size, 3 >(0) );
	parser_2_4.configure_rotameric_single_residue_dunbrack_library< 2, 4 >( lib_2_4, utility::fixedsizearray1< core::Size, 4 >(0) );
	parser_2_5.configure_rotameric_single_residue_dunbrack_library< 2, 5 >( lib_2_5, utility::fixedsizearray1< core::Size, 5 >(0) );
	parser_3_1.configure_rotameric_single_residue_dunbrack_library< 3, 1 >( lib_3_1, utility::fixedsizearray1< core::Size, 1 >(0) );
	parser_3_2.configure_rotameric_single_residue_dunbrack_library< 3, 2 >( lib_3_2, utility::fixedsizearray1< core::Size, 2 >(0) );
	parser_3_3.configure_rotameric_single_residue_dunbrack_library< 3, 3 >( lib_3_3, utility::fixedsizearray1< core::Size, 3 >(0) );
	parser_3_4.configure_rotameric_single_residue_dunbrack_library< 3, 4 >( lib_3_4, utility::fixedsizearray1< core::Size, 4 >(0) );
	parser_3_5.configure_rotameric_single_residue_dunbrack_library< 3, 5 >( lib_3_5, utility::fixedsizearray1< core::Size, 5 >(0) );
	parser_4_1.configure_rotameric_single_residue_dunbrack_library< 4, 1 >( lib_4_1, utility::fixedsizearray1< core::Size, 1 >(0) );
	parser_4_2.configure_rotameric_single_residue_dunbrack_library< 4, 2 >( lib_4_2, utility::fixedsizearray1< core::Size, 2 >(0) );
	parser_4_3.configure_rotameric_single_residue_dunbrack_library< 4, 3 >( lib_4_3, utility::fixedsizearray1< core::Size, 3 >(0) );
	parser_4_4.configure_rotameric_single_residue_dunbrack_library< 4, 4 >( lib_4_4, utility::fixedsizearray1< core::Size, 4 >(0) );
	parser_4_5.configure_rotameric_single_residue_dunbrack_library< 4, 5 >( lib_4_5, utility::fixedsizearray1< core::Size, 5 >(0) );
	parser_5_1.configure_rotameric_single_residue_dunbrack_library< 5, 1 >( lib_5_1, utility::fixedsizearray1< core::Size, 1 >(0) );
	parser_5_2.configure_rotameric_single_residue_dunbrack_library< 5, 2 >( lib_5_2, utility::fixedsizearray1< core::Size, 2 >(0) );
	parser_5_3.configure_rotameric_single_residue_dunbrack_library< 5, 3 >( lib_5_3, utility::fixedsizearray1< core::Size, 3 >(0) );
	parser_5_4.configure_rotameric_single_residue_dunbrack_library< 5, 4 >( lib_5_4, utility::fixedsizearray1< core::Size, 4 >(0) );
	parser_5_5.configure_rotameric_single_residue_dunbrack_library< 5, 5 >( lib_5_5, utility::fixedsizearray1< core::Size, 5 >(0) );
}

} //core
} //pack
} //dunbrack
