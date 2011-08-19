// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_pack_dunbrack_SingleResidueDunbrackLibrary_tmpl_hh
#define INCLUDED_core_pack_dunbrack_SingleResidueDunbrackLibrary_tmpl_hh

// Unit Headers
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

// Package Headers
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/chemical/AA.hh>

// C++ Headers
#include <algorithm>
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Residue.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.fwd.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>

// AUTO-REMOVED #include <fstream>

// Boost Headers
// AUTO-REMOVED #include <boost/cstdint.hpp>

namespace core {
namespace pack {
namespace dunbrack {

/*

/// DO NOT USE THIS CONSTRUCTOR
template < class T >
SingleResidueDunbrackLibraryConcrete< T >::SingleResidueDunbrackLibraryConcrete()
:
	SingleResidueDunbrackLibrary( chemical::aa_ala, 1, 1, false )
{}


template < class T >
SingleResidueDunbrackLibraryConcrete< T >::SingleResidueDunbrackLibraryConcrete(
	chemical::AA const aa_in,
	Size const nchi_aa_in,
	Size const nrots_per_bin,
	bool bClassicRotvectorAssignment
)
:
	SingleResidueDunbrackLibrary( aa_in, nchi_aa_in, nrots_per_bin, bClassicRotvectorAssignment ),
	rotset_( nrots_per_bin, nbins(), nbins() ), // column major ordering of FArrays
	max_rotprob_( nbins(), nbins(), (DunbrackReal) 0 )
{}

template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::set_chi_mean(
	Size phi_bin,
	Size psi_bin,
	Size which_rotamer,
	Size which_chi,
	DunbrackReal chi_mean_in )
{
	rotset_( which_rotamer, psi_bin, phi_bin ).chi_mean( which_chi, chi_mean_in );
}

template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::set_chi_sd(
	Size phi_bin,
	Size psi_bin,
	Size which_rotamer,
	Size which_chi,
	DunbrackReal chi_sd_in
)
{
	rotset_( which_rotamer, psi_bin, phi_bin ).chi_sd( which_chi, chi_sd_in );
}

template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::set_rotnum(
	Size phi_bin,
	Size psi_bin,
	Size which_rotamer,
	Size which_chi,
	Size rotnum_in )
{
	rotset_( which_rotamer, psi_bin, phi_bin ).rotnum( which_chi, rotnum_in );
}

template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::set_rotprob(
	Size phi_bin,
	Size psi_bin,
	Size which_rotamer,
	DunbrackReal rotprob_in
)
{
	//std::cout << "set rotprob: " << T::size << " rotno: ";
	//for ( Size ii = 1; ii <= T::size; ++ii ) {
	//	std::cout << rotset_( phi_bin, psi_bin, which_rotamer ).rotnum( ii ) << " ";
	//}
	//std::cout << std::endl;
	rotset_( which_rotamer, psi_bin, phi_bin ).rotamer_probability( rotprob_in );
	max_rotprob_( psi_bin, phi_bin ) = std::max( rotprob_in, max_rotprob_( psi_bin, phi_bin ) );
}

template < class T >
DunbrackReal
SingleResidueDunbrackLibraryConcrete< T >::get_max_rotprob(
	Size phi_bin,
	Size psi_bin
) const
{
	return max_rotprob_( psi_bin, phi_bin );
}

template < class T >
DunbrackRotamer< FOUR >
SingleResidueDunbrackLibraryConcrete< T >::find_rotamer(
	Size const phi_bin,
	Size const psi_bin,
	RotVector const & rot
) const
{
	Size correct_rot = 0;
	// replace this with templated functions
	for ( Size ii = 1; ii <= nrots_per_bin(); ++ii ) {
		bool all_match = true;
		for ( Size jj = 1; jj <= T::size; ++jj ) {
			if ( rotset_( ii, psi_bin, phi_bin ).rotnum( jj ) != (Size) rot[ jj ] ) {
				all_match = false;
				break;
			}
		}
		if ( all_match ) { correct_rot = ii; break; }
	}


	if ( correct_rot == 0 ) {

		std::cerr << "Failed to find rotamer: ";
		for ( Size ii = 1; ii <= T::size; ++ii ) {
			std::cerr << rot[ ii ] << " ";
		}
		std::cerr << std::endl << "Amongst options: ";
		for ( Size ii = 1; ii <= nrots_per_bin(); ++ii ) {
			for ( Size jj = 1; jj <= T::size; ++jj ) {
				std::cerr << rotset_( ii, psi_bin, phi_bin ).rotnum( jj ) << " ";
			}
			std::cerr << std::endl;
		}
		utility_exit();
	}

	return get_rotamer( phi_bin, psi_bin, correct_rot );
}

template < class T >
DunbrackRotamer< FOUR >
SingleResidueDunbrackLibraryConcrete< T >::find_rotamer(
	Size const phi_bin,
	Size const psi_bin,
	DunbrackRotamer< FOUR > const & dunrot
) const
{
	Size correct_rot = 0;
	// replace this with templated functions
	for ( Size ii = 1; ii <= nrots_per_bin(); ++ii ) {
		bool all_match = true;
		for ( Size jj = 1; jj <= T::size; ++jj ) {
			if ( rotset_( ii, psi_bin, phi_bin ).rotnum( jj ) != dunrot.rotnum( jj ) ) {
				all_match = false;
				break;
			}
		}
		if ( all_match ) { correct_rot = ii; break; }
	}

	if ( correct_rot == 0 ) {
		std::cerr << "Failed to find rotamer: ";
		for ( Size ii = 1; ii <= T::size; ++ii ) {
			std::cerr << dunrot.rotnum( ii ) << " ";
		}
		std::cerr << std::endl << "Amongst options: ";
		for ( Size ii = 1; ii <= nrots_per_bin(); ++ii ) {
			for ( Size jj = 1; jj <= T::size; ++jj ) {
				std::cerr << rotset_( ii, psi_bin, phi_bin ).rotnum( jj ) << " ";
			}
			std::cerr << std::endl;
		}
		utility_exit();
	}

	return get_rotamer( phi_bin, psi_bin, correct_rot );
}

template < class T >
DunbrackRotamer< FOUR >
SingleResidueDunbrackLibraryConcrete< T >::retrieve_rotamer(
	Size const phi_bin,
	Size const psi_bin,
	Size const which_rotamer
) const
{
	return get_rotamer( phi_bin, psi_bin, which_rotamer );
}


template < class T >
DunbrackRotamer< FOUR >
SingleResidueDunbrackLibraryConcrete< T >::get_rotamer(
	Size const phi_bin,
	Size const psi_bin,
	Size const which_rotamer
) const
{
	DunbrackRotamer< FOUR > rotamer;

	// replace this with templated functions
	for ( Size jj = 1; jj <= T::size; ++jj ) {
		rotamer.chi_mean( jj ) = rotset_( which_rotamer, psi_bin, phi_bin ).chi_mean( jj );
		rotamer.chi_sd( jj )   = rotset_( which_rotamer, psi_bin, phi_bin ).chi_sd(   jj );
		rotamer.rotnum( jj )   = rotset_( which_rotamer, psi_bin, phi_bin ).rotnum(   jj );
	}
	rotamer.rotamer_probability() = rotset_( which_rotamer, psi_bin, phi_bin ).rotamer_probability();
	return rotamer;
}

template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::write_to_binary(
	std::ofstream &out
) const
{
	using namespace boost;

	// 1. figure out how large each of the tables has to be
	// 2. allocate the c-style tables, fill them with data
	// 3. write the c-style tables in binary
	// 4. delete the c-style tables

	// 1.
	Size const nphipsibins = nbins() * nbins();
	Size nrotamers = nrots_per_bin();
	Size const rotamer_entries_needed = nphipsibins * nrotamers;
	Size const chi_entries_needed = rotamer_entries_needed * T::size;

	//SRDL_TR << "Writing " << nrotamers << " rotamers, with " << nchi_aa_ << " chi angles for " << aa_ << std::endl;

	// 2.
	DunbrackReal * prob_entries = new DunbrackReal[ rotamer_entries_needed ];
	boost::int32_t * rot_entries = new boost::int32_t[ chi_entries_needed ];
	DunbrackReal * chi_mean_entries = new DunbrackReal[ chi_entries_needed ];
	DunbrackReal * chi_sdev_entries = new DunbrackReal[ chi_entries_needed ];

	Size count_rots( 0 ), count_chis( 0 );
	for ( Size phi_bin = 1; phi_bin <= nbins(); phi_bin++ ) {
		for ( Size psi_bin = 1; psi_bin <= nbins(); psi_bin++ ) {
			for ( Size rot_id = 1; rot_id <= nrotamers; ++rot_id ) {
				DunbrackRotamer< T > const & rotamer = rotset_( rot_id, psi_bin, phi_bin );
				assert( count_rots < rotamer_entries_needed );
				prob_entries[ count_rots ] = rotamer.rotamer_probability();
				for ( Size chi_id = 1; chi_id <= T::size; ++chi_id ) {
					assert( count_chis < chi_entries_needed );
					rot_entries[ count_chis ]      = rotamer.rotnum( chi_id );
					chi_mean_entries[ count_chis ] = rotamer.chi_mean( chi_id );
					chi_sdev_entries[ count_chis ] = rotamer.chi_sd( chi_id );
					++count_chis;
				}
				++count_rots;
			}
		}
	}

	// 3.
	out.write( (char*) prob_entries, sizeof( DunbrackReal ) * rotamer_entries_needed );
	out.write( (char*) rot_entries,  sizeof( boost::int32_t ) * chi_entries_needed );
	out.write( (char*) chi_mean_entries,  sizeof( DunbrackReal ) * chi_entries_needed );
	out.write( (char*) chi_sdev_entries,  sizeof( DunbrackReal ) * chi_entries_needed );

	// 4.
	delete [] prob_entries;
	delete [] rot_entries;
	delete [] chi_mean_entries;
	delete [] chi_sdev_entries;

}


template < class T >
void
SingleResidueDunbrackLibraryConcrete< T >::read_from_binary(
	std::ifstream &in
)
{
	using namespace boost;

	// 1. figure out how large each of the tables has to be
	// 2. allocate the c-style tables,  read the binary file into the c-style tables
	// 3. store data from the c-style tables in the 3-D FArray
	// 4. delete the c-style tables

	// 1.
	Size nrotamers( nrots_per_bin() );

	Size const nphipsibins = nbins() * nbins();
	Size const rotamer_entries_needed = nphipsibins * nrotamers;
	Size const chi_entries_needed = rotamer_entries_needed * T::size;

	//SRDL_TR << "Reading " << nrotamers << " rotamers, with " << nchi_aa_ << " chi angles for " << aa_ << std::endl;


	// 2.
	DunbrackReal * prob_entries = new DunbrackReal[ rotamer_entries_needed ];
	boost::int32_t * rot_entries = new boost::int32_t[ chi_entries_needed ];
	DunbrackReal * chi_mean_entries = new DunbrackReal[ chi_entries_needed ];
	DunbrackReal * chi_sdev_entries = new DunbrackReal[ chi_entries_needed ];

	in.read( (char*) prob_entries, sizeof( DunbrackReal ) * rotamer_entries_needed );
	in.read( (char*) rot_entries,  sizeof( boost::int32_t ) * chi_entries_needed );
	in.read( (char*) chi_mean_entries,  sizeof( DunbrackReal ) * chi_entries_needed );
	in.read( (char*) chi_sdev_entries,  sizeof( DunbrackReal ) * chi_entries_needed );

	// 3.
	Size count_rots( 0 ), count_chis( 0 );
	for ( Size phi_bin = 1; phi_bin <= nbins(); phi_bin++ ) {
		for ( Size psi_bin = 1; psi_bin <= nbins(); psi_bin++ ) {
			for ( Size rot_id = 1; rot_id <= nrotamers; ++rot_id ) {
				assert( count_rots < rotamer_entries_needed );

				DunbrackRotamer< T > & rotamer( rotset_( rot_id, psi_bin, phi_bin ));
				//rotamer.rotamer_probability() = prob_entries[ count_rots ];
				set_rotprob(phi_bin, psi_bin, rot_id, prob_entries[ count_rots ]);
				for ( Size chi_id = 1; chi_id <= T::size; ++chi_id ) {
					assert( count_chis < chi_entries_needed );
					rotamer.rotnum( chi_id )   = rot_entries[ count_chis ];
					rotamer.chi_mean( chi_id ) = chi_mean_entries[ count_chis ];
					rotamer.chi_sd( chi_id ) = chi_sdev_entries[ count_chis ];
					++count_chis;
				}
				++count_rots;
			}
		}
	}

	// 4.
	delete [] prob_entries;
	delete [] rot_entries;
	delete [] chi_mean_entries;
	delete [] chi_sdev_entries;


}

*/

} // namespace dunbrack
} // namespace scoring
} // namespace core


#endif
