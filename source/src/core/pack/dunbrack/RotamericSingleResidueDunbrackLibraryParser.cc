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
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

//Numeric includes
#include <numeric/angle.functions.hh>

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
	chi_std_devs_(),
	is_canonical_dun02_library_(false)
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
	chi_std_devs_(src.chi_std_devs_),
	is_canonical_dun02_library_(src.is_canonical_dun02_library_)
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
	is_canonical_dun02_library_ = false;

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

			// Skip checks and corrections for the old Dun02 libraries for canonical amino acids.
			// Note: When we finally remove score12 support, this should also be removed.
			if ( core::chemical::is_canonical_L_aa_or_gly( aa ) && !is_shapovalov_file ) is_canonical_dun02_library_ = true;
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
				do_all_checks_and_corrections( infile.filename() );
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

	do_all_checks_and_corrections( infile.filename() );

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

/// @brief Determine the remapping for the indices of the rotamer wells read from disk to ensure that they
/// have the order consistent with the Rosetta convention.
/// @details Rotamers should be in order from lowest to highest in the range [0, 360).  Returns true if this is NOT the case, false otherwise.
/// If true, the rotamer_well_reordering map (which is cleared by this operation) is populated with a remapping that sorts the rotamers properly.
/// @note This uses the first backbone bin to determine the reordering.  After applying this to all rotamer wells, it is necessary
/// to check that no further reordering is needed.  The filename parameter is only used for error messages.
bool
RotamericSingleResidueDunbrackLibraryParser::determine_rotamer_well_order(
	utility::vector1< std::map< core::Size, core::Size > > & rotamer_well_reordering,
	std::string const &filename
) const {
	core::Size const n_rotwells( rotwells_.size() );
	runtime_assert_string_msg( n_rotwells > 0, "Error in RotamericSingleResidueDunbrackLibraryParser::determine_rotamer_well_order(): No rotamer wells have yet been stored!" );

	bool reordering_needed( false );

	utility::vector1< core::Size > first_bb_bin_indices; //List of indices of rotwells_ entries that correspond to the first backbone bin.
	utility::vector1< core::Size > max_rotwell_indices; //The maximum values of the rotamer well indices for each chi.
	max_rotwell_indices.resize( rotwells_[1].size(), 0 );
	for ( core::Size i(1); i<=n_rotwells; ++i ) {
		if ( backbone_torsions_[i] == backbone_torsions_[1] ) first_bb_bin_indices.push_back(i);
		for ( core::Size j(1); j<=rotwells_[i].size(); ++j ) {
			debug_assert(rotwells_[i].size() == max_rotwell_indices.size());
			if ( max_rotwell_indices[j] < rotwells_[i][j] ) max_rotwell_indices[j] = rotwells_[i][j]; //Find maximum rotamer well indices
		}
	}
	core::Size const n_rotwells_in_first_bb_bin( first_bb_bin_indices.size() );

	rotamer_well_reordering.clear();
	rotamer_well_reordering.resize( rotwells_[1].size() );

	// For each rotamer well for each chi in the first backbone bin:
	// (a) Find the minimum of the chimean values.
	// (b) Sort the rotamer well indices by the min chimean values, in the range [0, 360).
	// (c) Set up the mapping.
	for ( core::Size ichi(1), ichimax( rotwells_[1].size() ); ichi<=ichimax; ++ichi ) {
		// (a) Find the minimum of the chimean values:
		std::map< core::Size, core::Real > chimean_mins;
		core::Size n_rotwells(0);
		for ( core::Size i(1); i<=n_rotwells_in_first_bb_bin; ++i ) {
			core::Size const cur_rotwell(rotwells_[i][ichi]);
			if ( cur_rotwell == 0 ) continue; //zero rotwell -- indicates rotwells not defined for this chi.
			if ( !chimean_mins.count( cur_rotwell ) ) ++n_rotwells;
			if ( !chimean_mins.count( cur_rotwell ) || ( chimean_mins.at(cur_rotwell) > numeric::nonnegative_principal_angle_degrees( chimeans_[i][ichi] ) )  ) {
				chimean_mins[ cur_rotwell ] = numeric::nonnegative_principal_angle_degrees( chimeans_[i][ichi] );
			}
		}

		if ( n_rotwells == 0 ) continue; // No rotamer wells for this chi.

		core::Real biggestval(0); bool first(true);
		for ( std::map< core::Size, core::Real >::iterator mapit( chimean_mins.begin() ); mapit != chimean_mins.end(); ++mapit ) {
			if ( first || mapit->second > biggestval ) {
				first = false;
				biggestval = mapit->second; //Figure out the biggest value in the list.
			}
		}

		// (b) Sort the rotamer well indices by the min chimean values.  We'll use an inefficent selection sort, since the list is guaranteed to be small:
		for ( core::Size i(1); i<=n_rotwells; ++i ) { //Outer loop for n^2 selection sort.
			core::Real minval(0); core::Size minj(0); first=true /*Reuse this var from above*/;
			for ( core::Size j(1); j<=n_rotwells; ++j ) { //Inner loop for n^2 selection sort.
				runtime_assert_string_msg( chimean_mins.count(j), "Error in RotamericSingleResidueDunbrackLibraryParser::determine_rotamer_well_order(): Rotamer file " + filename +  " does not have continuously-numbered rotamer wells!" ); //Should always be true, unless file is badly formatted.
				if ( !value_is_in_map( rotamer_well_reordering[ichi], j ) /*Not already sorted*/ && (first || chimean_mins.at(j) < minval) ) { // Find the smallest value that isn't already sorted.
					minval = chimean_mins.at(j);
					minj = j;
					first=false;
				}
			} //Inner loop for n^2 selection sort.

			// (c) Set up the mapping:
			debug_assert( !rotamer_well_reordering[ichi].count(i) );
			rotamer_well_reordering[ichi][i] = minj;
			if ( minj != i ) reordering_needed = true;
		} //Outer loop for n^2 selection sort.
	}

	if ( TR.Debug.visible() ) {
		for ( core::Size i(1), imax(rotamer_well_reordering.size()); i<=imax; ++i ) {
			TR.Debug << "CHI" << i << ":";
			for ( core::Size j(1), jmax(rotamer_well_reordering[i].size()); j<=jmax; ++j ) {
				TR.Debug << "\t" << rotamer_well_reordering[i].at(j);
			}
			TR.Debug << std::endl;
		}
		TR.Debug.flush();
	}

	return reordering_needed;
}

/// @brief This does a simple check of all of the rotamer well order, to determine whether it's consistent with
/// the Rosetta convention.
/// @details Rotamers should be in order from lowest to highest in the range [0, 360).  Returns true if this is NOT the case, false otherwise.
/// @note This checks ALL backbone bins.
bool
RotamericSingleResidueDunbrackLibraryParser::check_rotamer_well_order() const {

	runtime_assert_string_msg( rotwells_.size() > 0, "Error in RotamericSingleResidueDunbrackLibraryParser::check_rotamer_well_order(): No rotamer wells have yet been stored!" );

	utility::vector1< utility::vector1< core::Real > > backbone_torsions_completed;
	utility::vector1< core::Size > n_rotwells_for_each_chi(rotwells_[1].size(), 0);
	bool first_bb_bin(true);

	for ( core::Size ii(1), iimax(rotwells_.size()); ii<=iimax; ++ii ) {
		core::Size const n_chi( rotwells_[ii].size() );

		bool continue_on(false);
		if ( !backbone_torsions_completed.empty() ) {
			for ( core::Size i(1), imax(backbone_torsions_completed.size()); i<=imax; ++i ) {
				if ( backbone_torsions_[ii] == backbone_torsions_completed[i] ) { //We've already done this backbone torsion well.
					continue_on = true;
					break;
				}
			}
		}
		if ( continue_on ) continue;
		backbone_torsions_completed.push_back(backbone_torsions_[ii]);

		utility::vector1< core::Size > bb_bin_indices; //List of indices of rotwells_ entries that correspond to the current backbone bin.
		for ( core::Size i(1), imax(rotwells_.size()); i<=imax; ++i ) {
			if ( backbone_torsions_[i] == backbone_torsions_[ii] ) {
				bb_bin_indices.push_back(i);
				if ( first_bb_bin ) {
					for ( core::Size j(1); j<=n_chi; ++j ) {
						if ( rotwells_[i][j] > n_rotwells_for_each_chi[j] ) n_rotwells_for_each_chi[j] = rotwells_[i][j];
					}
				}
			}
		}
		first_bb_bin = false;
		core::Size const n_entries_for_this_bin( bb_bin_indices.size() );

		//For each rotamer in the Nth backbone bin, check whether every lower-index rotamer has lower mean chi values:
		for ( core::Size i(1); i<=n_entries_for_this_bin; ++i ) { //Loop over all rotamer wells in this backbone bin
			core::Size const iprime( bb_bin_indices[i] );
			for ( core::Size j(1); j<=n_entries_for_this_bin; ++j ) { //Inner loop over all rotamer wells in this backbone bin, for comparison
				core::Size const jprime( bb_bin_indices[j] );
				if ( i==j ) continue; //No self-comparisons
				for ( core::Size k(1), kmax(rotwells_[iprime].size()); k<=kmax; ++k ) { //Loop over all chis
					if ( rotwells_[iprime][k] == 1 ) continue;
					if ( rotwells_[iprime][k] > rotwells_[jprime][k] ) { //If we're comparing to a lower-index rotamer, given the remapping established already...
						if ( numeric::nonnegative_principal_angle_degrees( chimeans_[iprime][k] ) < numeric::nonnegative_principal_angle_degrees( chimeans_[jprime][k] ) ) { //Then the chimeans of the upper rotamer should be greater than those of the lower.
							if ( rotwells_[iprime][k] != n_rotwells_for_each_chi[k] && rotwells_[jprime][k] !=1 ) return true; //The highest or lowest rotamer wells should be allowed to wrap around.
						}
					}
				}
			}
		}

	}

	return false;
}

/// @brief Given a map defining a reordering of rotamer wells, reorder rotamer wells.
void
RotamericSingleResidueDunbrackLibraryParser::update_rotamer_well_order(
	utility::vector1< std::map< core::Size, core::Size > > const & rotamer_well_reordering
) {
	for ( core::Size i(1), imax( rotwells_.size() ); i<=imax; ++i ) { //Loop through all rotamer wells
		utility::vector1< core::Size > & cur_rotwells( rotwells_[i] );
		debug_assert( cur_rotwells.size() == rotamer_well_reordering.size() );
		for ( core::Size j(1), jmax( cur_rotwells.size() ); j<=jmax; ++j ) { //Loop through all chis
			if ( rotamer_well_reordering[j].count( cur_rotwells[j] ) ) {
				cur_rotwells[j] = rotamer_well_reordering[j].at( cur_rotwells[j] );
			}
		}
	}
}

/// @brief Performs a number of checks and corrections.
/// @details Calls check_correct_vector_lengths(), determine_rotamer_well_order(), check_rotamer_well_order(), and update_rotamer_well_order(), and , currently.  Additional checks
/// and corrections may be added in the future.
/// @param[in] filename The name of the rotamer file currently being read.  This is only used for error messages.
void
RotamericSingleResidueDunbrackLibraryParser::do_all_checks_and_corrections( std::string const & filename ) {
	check_correct_vector_lengths();

	//Checks for rotamer well order:
	if ( !is_canonical_dun02_library_ ) {
		utility::vector1< std::map <core::Size, core::Size> > rotamer_well_reordering;
		bool const rotamer_well_order_bad( determine_rotamer_well_order( rotamer_well_reordering, filename ) );
		if ( rotamer_well_order_bad ) {
			TR.Warning << TR.Red << "Rotamer file " << filename << " has incorrect rotamer well assignments.  The Rosetta convention is to number rotamer wells from lowest to highest in the range [0,360).  Attempting to correct automatically." << TR.Reset << std::endl;
			update_rotamer_well_order( rotamer_well_reordering );
			if ( check_rotamer_well_order() ) {
				TR.Warning << TR.Red << "Rotamer file " << filename << " has inconsistent rotamer well assignments even after reordering.  This may or may not indicate a badly formatted file." << TR.Reset << std::endl;
			}
		}
	} else {
		bool const old_dun02_library( is_old_canonical_dun02_library( filename ) );
		bool const beta_nov16_library( is_canonical_beta_nov16_library(filename) );
		bool const talaris_library( is_canonical_talaris_library(filename) );
		if ( !old_dun02_library && !beta_nov16_library && !talaris_library ) {
			TR.Warning << TR.Red << "Rotamer file " << filename  << " appears to be a non-canonical rotamer library, of Dunbrack2002 format, for a canonical amino acid.  This can happen, for example, when a canonical amino acid is modified (e.g. N-methylated).  It will be assumed that rotamer well identities match those for the corresponding canonical amino acid.  Skipping rotamer order corrections."  << TR.Reset << std::endl;
		} else if ( (beta_nov16_library || talaris_library) && !old_dun02_library ) {
			TR.Debug << "Rotamer file " << filename << " is a canonical rotamer file for ";
			if ( beta_nov16_library ) {
				TR.Debug << "the new beta_nov16 ";
			} else {
				TR.Debug << "the old talaris ";
			}
			TR.Debug << "energy function (which uses Dunbrack2002-formatted rotamer files).  Skipping rotamer well order corrections." << std::endl;
			if ( talaris_library ) {
				TR.Debug << "Please consider switching to a newer energy function, since support for the talaris energy function will not be maintained indefinitely." << std::endl;
			}
		} else {
			TR.Debug << "Rotamer file " << filename << " is a canonical rotamer file for the old Dunbrack2002 libraries.  Skipping rotamer well order corrections -- but please note that these libraries will soon be deprecated.  Please consider switching to the current default scoring function." << std::endl;
		}
	}
}

/// @brief Given a filename, return true if this is a talaris library for a canonical amino acid, false otherwise.
bool
RotamericSingleResidueDunbrackLibraryParser::is_canonical_talaris_library(
	std::string const &filename
) const {
	std::string const path( utility::pathname( filename ) );
	if ( path.find( "/rotamer/ExtendedOpt1-5/" ) == std::string::npos ) return false;

	std::string const filename_stripped( utility::filename( filename ) );
	return
		!filename_stripped.compare( "cys.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "ile.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "lys.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "leu.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "met.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "pro.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "arg.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "ser.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "thr.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "val.bbdep.rotamers.lib.gz" );
}


/// @brief Given a filename, return true if this is a beta_nov16 library for a canonical amino acid, false otherwise.
bool
RotamericSingleResidueDunbrackLibraryParser::is_canonical_beta_nov16_library(
	std::string const &filename
) const {
	std::string const path( utility::pathname( filename ) );
	if ( path.find( "/rotamer/beta_nov2016/" ) == std::string::npos ) return false;

	std::string const filename_stripped( utility::filename( filename ) );
	return
		!filename_stripped.compare( "cys.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "ile.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "lys.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "leu.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "met.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "pro.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "arg.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "ser.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "thr.bbdep.rotamers.lib.gz" ) ||
		!filename_stripped.compare( "val.bbdep.rotamers.lib.gz" );
}

/// @brief Given a filename, return true if this is an old Dunbrack 2002 library for a canonical amino acid, false otherwise.
bool
RotamericSingleResidueDunbrackLibraryParser::is_old_canonical_dun02_library(
	std::string const &filename
) const {
	return !( utility::filename( filename ).compare("bbdep02.May.sortlib-correct.12.2010") ) || !( utility::filename( filename ).compare("bbdep02.May.sortlib") );
}

/// @brief Given a Size->Size map, determine whether there exists a key that maps to a given value.
bool
RotamericSingleResidueDunbrackLibraryParser::value_is_in_map(
	std::map< core::Size, core::Size > const &the_map,
	core::Size the_value
) const {
	for ( std::map<core::Size, core::Size>::const_iterator it( the_map.begin() ); it!=the_map.end(); ++it ) {
		if ( it->second == the_value ) return true;
	}
	return false;
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
