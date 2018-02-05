// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.hh
/// @brief A helper class to assist in parsing rotamer libraries.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_hh

#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.fwd.hh>

// Core headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

namespace core {
namespace pack {
namespace dunbrack {

/// @brief A helper class to assist in parsing rotamer libraries.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
class RotamericSingleResidueDunbrackLibraryParser : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor.
	/// @param[in] num_mainchain_torsions The number of mainchain torsions on which this rotamer library depends.
	/// @param[in] max_possible_chis The maximum number of chis allowed in this rotamer library.  This is the number of chi columns in the file (usually 4).
	/// @param[in] max_possible_rotamers The maximum number of rotamers that this rotamer library can define.
	RotamericSingleResidueDunbrackLibraryParser( core::Size const num_mainchain_torsions, core::Size const max_possible_chis, core::Size const max_possible_rotamers );

	RotamericSingleResidueDunbrackLibraryParser(RotamericSingleResidueDunbrackLibraryParser const & src);

	virtual ~RotamericSingleResidueDunbrackLibraryParser();

	RotamericSingleResidueDunbrackLibraryParserOP
	clone() const;

	/// @brief Clear all of the data read from a file.
	/// @details Also resets read_file_was_called_ to false.
	void clear_imporated_data();

	/// @brief Read an input stream to set up this object.
	/// @details The RotamericSingleResidueDunbrackLibraryParser stores the information from the rotamer library in a very one-to-one
	/// format, and can subsequently be used to configure a RotamericSingleResidueDunbrackLibrary object (requiring some interpretation
	/// of the data).
	/// @note We might be starting the read in the middle of the file, and the three-letter code from the first line might already have been read.
	/// @returns The name of the next amino acid specified in the input stream, or alternatively an empty string if EOF is reached.
	std::string read_file( utility::io::izstream & infile, bool first_line_three_letter_code_already_read, core::chemical::AA const aa );

	/// @brief Based on the material read from the file and stored in this object, configure
	/// a dunbrack library.
	/// @details read_file() must be called first.
	/// @param[in] library Pass in *this from a RotamericSingleResidueDunbrackLibrary member function.
	/// @param[in] n_bb_bins Pass in N_BB_BINS.
	/// @note Implemented in RotamericSingleResidueDunbrackLibraryParser.tmpl.hh.
	template< core::Size T, core::Size N >
	void
	configure_rotameric_single_residue_dunbrack_library(
		RotamericSingleResidueDunbrackLibrary< T, N > &library,
		utility::fixedsizearray1< Size, N > const & n_bb_bins
	) const;

	/// @brief Calculate the product of elements in a fixedsizearray.
	template< Size N >
	inline Size calc_product( utility::fixedsizearray1< Size, N > factors ) const {
		Size x = factors[ 1 ];
		for ( Size i = 2; i <= N; ++i ) x *= factors[ i ];
		return x;
	}


private:
	/////PRIVATE MEMBER FUNCTIONS/////

	/// @brief When the first non-comment line is read from a rotamer file, check whether there's an extra column (signifying that this is a
	/// Shapovalov file, which has an extra column in the middle.)
	/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
	bool check_for_extra_column( utility::io::izstream & infile, bool const first_line_three_letter_code_already_read ) const;

	/// @brief The legacy reader had a check for mainchain torsions equal to -180, in which case the last rotamer read was discarded.
	/// Maintaining that here.
	/// @returns Returns true for correction applied (i.e. last rotamer deleted), false otherwise.
	bool correct_for_duplicated_mainchain_data();

	/// @brief Delete the last entry in all of the data vectors.
	void remove_last_entry_in_vectors();

	/// @brief Confirm that all data vectors have the same length.  Throw an error if they don't.
	void check_correct_vector_lengths() const;

	/// @brief Determine the remapping for the indices of the rotamer wells read from disk to ensure that they
	/// have the order consistent with the Rosetta convention.
	/// @details Rotamers should be in order from lowest to highest in the range [0, 360).  Returns true if this is NOT the case, false otherwise.
	/// If true, the rotamer_well_reordering map (which is cleared by this operation) is populated with a remapping that sorts the rotamers properly.
	/// @note This uses the first backbone bin to determine the reordering.  After applying this to all rotamer wells, it is necessary
	/// to check that no further reordering is needed.  The filename parameter is only used for error messages.
	bool determine_rotamer_well_order( utility::vector1< std::map< core::Size, core::Size > > & rotamer_well_reordering, std::string const &filename ) const;

	/// @brief This does a simple check of all of the rotamer well order, to determine whether it's consistent with
	/// the Rosetta convention.
	/// @details Rotamers should be in order from lowest to highest in the range [0, 360).  Returns true if this is NOT the case, false otherwise.
	/// @note This checks ALL backbone bins.
	bool check_rotamer_well_order() const;

	/// @brief Given a map defining a reordering of rotamer wells, reorder rotamer wells.
	void update_rotamer_well_order( utility::vector1< std::map< core::Size, core::Size > > const & rotamer_well_reordering );

	/// @brief Performs a number of checks and corrections.
	/// @details Calls check_correct_vector_lengths(), determine_rotamer_well_order(), check_rotamer_well_order(), and update_rotamer_well_order(), and , currently.  Additional checks
	/// and corrections may be added in the future.
	/// @param[in] filename The name of the rotamer file currently being read.  This is only used for error messages.
	void do_all_checks_and_corrections( std::string const & filename );

	/// @brief Given a filename, return true if this is a talaris library for a canonical amino acid, false otherwise.
	bool is_canonical_talaris_library( std::string const &filename) const;

	/// @brief Given a filename, return true if this is a beta_nov16 library for a canonical amino acid, false otherwise.
	bool is_canonical_beta_nov16_library( std::string const &filename) const;

	/// @brief Given a filename, return true if this is an old Dunbrack 2002 library for a canonical amino acid, false otherwise.
	bool is_old_canonical_dun02_library( std::string const &filename ) const;

	/// @brief Given a Size->Size map, determine whether there exists a key that maps to a given value.
	bool value_is_in_map( std::map< core::Size, core::Size > const &the_map, core::Size the_value ) const;

private:

	/// @brief This function forces the instantiation of virtual templated methods in the derived classes.
	/// Functions like this one are necessary when combining polymorphism and templates.  Though
	/// these functions must be compiled, they need never be called. Do not call this function.
	/// @details THIS FUNCTION SHOULD NEVER BE CALLED!
	void hacky_parser_init_workaround();

private:
	/////PRIVATE DATA/////

	/// @brief Was the read_file() method called?
	bool read_file_was_called_;

	/// @brief Number of mainchain torsions upon which this rotamer library depends.
	core::Size num_mainchain_torsions_;

	/// @brief Max number of chis allowed.
	core::Size max_possible_chis_;

	/// @brief The maximum number of possible rotamers.
	core::Size max_possible_rotamers_;

	/////DATA READ FROM FILE BELOW THIS POINT//////

	/// @brief Backbone torsions for each rotamer.
	utility::vector1< utility::vector1< core::Real > > backbone_torsions_;

	/// @brief Counts for each rotamer.  Not actually used.
	utility::vector1< core::Size > counts_;

	/// @brief Rotamer indices for each rotamer.  Only used for warnings.
	utility::vector1< utility::vector1 < core::Size > > rotwells_;

	/// @brief Probabilities for each rotamer.  Used in scoring.
	utility::vector1< core::Real > probabilities_;

	/// @brief The chimean values, which define the centres of each rotamer.
	utility::vector1< utility::vector1 < core::Real > > chimeans_;

	/// @brief The chi standard deviation values, which define the ramping of each chi well.
	utility::vector1< utility::vector1 < core::Real > > chi_std_devs_;

	/// @brief Is this library a canonical amino acid's Dunbrack2002 library?
	/// @details If it is, we need to skip the correction step, since the canonical Dun02 library has
	/// wonky rotamer well orders that are corrected some other way.
	bool is_canonical_dun02_library_;

};

} //core
} //pack
} //dunbrack



#endif //INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_hh
