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
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
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
	template< core::Size T, core::Size N >
	void
	configure_rotameric_single_residue_dunbrack_library(
		RotamericSingleResidueDunbrackLibrary< T, N > &library,
		utility::fixedsizearray1< Size, N > const & n_bb_bins
	) const {
		runtime_assert_string_msg( read_file_was_called_, "Error in core::pack::dunbrack::RotamericSingleResidueDunbrackLibraryParser::configure_rotameric_single_residue_dunbrack_library(): The read_file() method must be called first!" );
		debug_assert( backbone_torsions_.size() == rotwells_.size() ); //Should be true.

		core::Size const num_bb_bins( calc_product( n_bb_bins ) );

		typename utility::vector1< DunbrackRotamer< T, N > > first_bb_bin_data;
		first_bb_bin_data.reserve( library.n_possible_rots() );

		bool finished_first_bb_bin(false); //Have we finished the first backbone bin?
		bool very_first_rotamer( true ); //Is this the very first rotamer that we're parsing?

		utility::fixedsizearray1< core::Size, N > last_bb_bin(0);
		utility::vector1< core::Size > bb_bin( N, 0 );
		core::Size count_in_this_bb_bin( 1 );

		for ( core::Size i(1), imax(backbone_torsions_.size()); i<=imax; ++i ) { //Loop through all rotamers
			utility::vector1< core::Real > const & bb_torsions( backbone_torsions_[i] );

			library.get_bb_bins( bb_torsions, bb_bin );

			if ( finished_first_bb_bin ) {
				Size const packed_rotno( library.rotwell_2_packed_rotno( rotwells_[i] ) );

				if ( packed_rotno < 1 || packed_rotno > library.n_packed_rots() ) {
					std::cerr << "ERROR in converting rotwell to packed rotno: ";
					for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwells_[i][ ii ];
					std::cerr << " " << packed_rotno << std::endl;
					utility_exit();
				}

				PackedDunbrackRotamer< T, N > rotamer( chimeans_[i], chi_std_devs_[i], probabilities_[i], packed_rotno );

				for ( Size ii = 1; ii <= N; ++ii ) {
					if ( bb_bin[ ii ] != last_bb_bin[ ii ] ) {
						count_in_this_bb_bin = 1; break;
					}
				}
				if ( count_in_this_bb_bin > library.n_packed_rots() ) break;

				core::Size bb_bin_index = make_index< N >( n_bb_bins, bb_bin );

				library.rotamers()( bb_bin_index, count_in_this_bb_bin ) = rotamer;

				library.packed_rotno_2_sorted_rotno()( bb_bin_index, packed_rotno ) = count_in_this_bb_bin;

				++count_in_this_bb_bin;
				for ( Size ii = 1; ii <= N; ++ii ) last_bb_bin[ ii ] = bb_bin[ ii ];
			} else {
				bool notdone = false;
				for ( Size ii = 1; ii <= N; ++ii ) {
					if ( bb_bin[ ii ] != last_bb_bin[ ii ] ) {
						notdone = true; break;
					}
				}

				if ( !very_first_rotamer && notdone ) {
					// We have now read all rotamers from this phi/psi bin.
					// 1. Declare all existing rotwells encountered
					// 2. Allocate space for rotamer data
					// 3. Store data from first rotwell.
					// 4. Store the data from the rotamer we just read.

					// 1.
					library.declare_all_existing_rotwells_encountered();
					// 2.
					library.rotamers().dimension( num_bb_bins, library.n_packed_rots() );//, PackedDunbrackRotamer< T, N >() );
					library.packed_rotno_2_sorted_rotno().dimension( num_bb_bins, library.n_packed_rots() );
					// 3.
					utility::vector1< Size > first_bin_rotwell( DUNBRACK_MAX_SCTOR, 0 );
					for ( Size ii = 1; ii <= first_bb_bin_data.size(); ++ii ) {
						for ( Size jj = 1; jj <= T; ++jj ) first_bin_rotwell[ jj ] = first_bb_bin_data[ ii ].rotwell( jj );
						Size const packed_rotno = library.rotwell_2_packed_rotno( first_bin_rotwell );

						if ( packed_rotno < 1 || packed_rotno > library.n_packed_rots() ) {
							std::cerr << "ERROR in converting rotwell to packed rotno: ";
							for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwells_[i][ ii ];
							std::cerr << " " << packed_rotno << std::endl;
							utility_exit();
						}

						Size bb_bin_index = make_index< N >( n_bb_bins, last_bb_bin );
						library.rotamers()( bb_bin_index, ii ) = PackedDunbrackRotamer< T, N >( first_bb_bin_data[ ii ], packed_rotno );
						library.rotamers()( bb_bin_index, ii ).rotamer_probability() = first_bb_bin_data[ ii ].rotamer_probability();
						library.packed_rotno_2_sorted_rotno()( bb_bin_index, packed_rotno ) = ii;
					}

					// 4.
					debug_assert( count_in_this_bb_bin == 1 );
					Size const packed_rotno( library.rotwell_2_packed_rotno( rotwells_[i] ) );
					if ( packed_rotno < 1 || packed_rotno > library.n_packed_rots() ) {
						std::cerr << "ERROR in converting rotwell to packed rotno: ";
						for ( Size ii = 1; ii <= T; ++ii ) std::cerr << " " << rotwells_[i][ ii ];
						std::cerr << " " << packed_rotno << std::endl;
						utility_exit();
					}
					PackedDunbrackRotamer< T, N > rotamer( chimeans_[i], chi_std_devs_[i], probabilities_[i], packed_rotno );
					Size bb_bin_index = make_index< N >( n_bb_bins, bb_bin );

					library.rotamers()( bb_bin_index, count_in_this_bb_bin ) = rotamer;
					library.packed_rotno_2_sorted_rotno()( bb_bin_index, packed_rotno ) = count_in_this_bb_bin;
					++count_in_this_bb_bin;
					for ( Size bbi = 1; bbi <= N; ++bbi ) last_bb_bin[ bbi ] = bb_bin[ bbi ];

					finished_first_bb_bin = true;
				} else {
					very_first_rotamer = false;
					library.mark_rotwell_exists( rotwells_[i] );
					first_bb_bin_data.push_back( DunbrackRotamer< T, N >( chimeans_[i], chi_std_devs_[i], probabilities_[i], rotwells_[i] ) );
					first_bb_bin_data[first_bb_bin_data.size()].rotamer_probability() = probabilities_[i];
					for ( Size bbi = 1; bbi <= N; ++bbi ) last_bb_bin[ bbi ] = bb_bin[ bbi ];
				}
			}

		}
	}

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
};

} //core
} //pack
} //dunbrack



#endif //INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_hh
