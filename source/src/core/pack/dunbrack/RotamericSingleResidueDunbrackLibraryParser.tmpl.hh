// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.tmpl.hh
/// @brief A helper class to assist in parsing rotamer libraries.  Template functions are declared here.
/// @author Vikram K. Mulligan (vmullig@uw.edu).


#ifndef INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_tmpl_hh
#define INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_tmpl_hh

#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.fwd.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.hh>
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

/// @brief Based on the material read from the file and stored in this object, configure
/// a dunbrack library.
/// @details read_file() must be called first.
/// @param[in] library Pass in *this from a RotamericSingleResidueDunbrackLibrary member function.
/// @param[in] n_bb_bins Pass in N_BB_BINS.
/// @note Implemented in RotamericSingleResidueDunbrackLibraryParser.tmpl.hh.
template< core::Size T, core::Size N >
void
RotamericSingleResidueDunbrackLibraryParser::configure_rotameric_single_residue_dunbrack_library(
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

} //core
} //pack
} //dunbrack



#endif //INCLUDED_core_pack_dunbrack_RotamericSingleResidueDunbrackLibraryParser_tmpl_hh
