// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/AtomVDW.hh
/// @brief
/// @author Rhiju Das


// Unit Headers
#include <core/scoring/rna/RNA_AtomVDW.hh>

// Package headers

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.rna.RNA_AtomVDW" );

namespace core {
namespace scoring {
namespace rna {

Size
rna_residue_name_to_num( char const c )
{
	if ( c == 'a' ) return 1;
	if ( c == 'c' ) return 2;
	if ( c == 'g' ) return 3;
	if ( c == 'u' ) return 4;
	if ( c == 't' ) return 4;
	if ( c == 'Z' ) return 5; // Mg(2+)
	tr << "What is this? " << c << std::endl;
	utility_exit_with_message( "Asked for rna_residue_name_to_num for unknown residue_name" );
	return 0;
}

//Why doesn't this helper function already exist in vector class?
template
< class T >
Size
get_position_in_vector( utility::vector1< T > & vec, T const & element )
{
	Size count( 1 );
	for ( auto const & elem : vec ) {
		if ( elem == element ) return count;
		count++;
	}
	vec.push_back( element );
	return count;
}


/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
RNA_AtomVDW::RNA_AtomVDW()
{
	//Read in data file, and fill in private data.
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/rna/rna_atom_vdw.txt" );

	rna_vdw_parameter_.dimension( 13, 13, 5, 5 );  // 5 = a,c,g,u,Z (Mg2+)
	rna_vdw_parameter_ = 0.0; // zero everything out.

	if ( !stream.good() ) utility_exit_with_message( "Unable to open rna_scoring/AtomVDW/atom_vdw.txt!" );

	// read the entire file and figure out what atom_types are present and in what order
	utility::vector1< std::string > lines;

	std::string line;
	std::string atom_name1, atom_name2;
	char which_residue1, which_residue2;
	Real input_bump_parameter;
	while ( getline( stream, line ) ) {
		lines.push_back( line );
		std::istringstream l( line );
		l >> which_residue1 >> which_residue2 >> atom_name1 >> atom_name2 >> input_bump_parameter;

		Size const pos1 = get_position_in_vector( rna_vdw_atom_[which_residue1], atom_name1 );
		Size const pos2 = get_position_in_vector( rna_vdw_atom_[which_residue2], atom_name2 );

		runtime_assert( pos1 <= rna_vdw_parameter_.size1() );
		runtime_assert( pos1 <= rna_vdw_parameter_.size2() );

		rna_vdw_parameter_( pos1,
			pos2,
			rna_residue_name_to_num( which_residue1 ),
			rna_residue_name_to_num( which_residue2 ) ) = input_bump_parameter;
		//perhaps we should explicitly force symmetry here?
	}

	// for ( Size i = 1; i <= 9; i++ ) {
	//  for ( Size j = 1; j <= 9; j++ ) {
	//   for ( Size m = 1; m <= 5; m++ ) {
	//    for ( Size n = 1; n <= 5; n++ ) {
	//     std::cout << "BUMP PARAM: " << i << " " << j << " " << m << " " << n << " " << rna_vdw_parameter_( i, j, m, n ) << std::endl;
	//    }
	//   }
	//  }
	// }
}

//////////////////////////////////////////////
utility::vector1 < std::string > const
RNA_AtomVDW::vdw_atom_list( chemical::ResidueType const & rt ) const
{
	using namespace core::chemical;
	char which_nucleotide = aa_unp;
	if ( rt.aa() == na_rad || rt.na_analogue() == na_rad ) {
		which_nucleotide = 'a';
	} else if ( rt.aa() == na_rcy || rt.na_analogue() == na_rcy ) {
		which_nucleotide = 'c';
	} else if ( rt.aa() == na_rgu || rt.na_analogue() == na_rgu ) {
		which_nucleotide = 'g';
	} else if ( rt.aa() == na_ura || rt.na_analogue() == na_ura ) {
		which_nucleotide = 'u';
	} else if ( rt.name1() == 't' ) {
		which_nucleotide = 't';
	}
	
	auto iter = rna_vdw_atom_.find( which_nucleotide );

	if ( iter == rna_vdw_atom_.end() ) {
		if ( which_nucleotide == 't' ) {
			iter = rna_vdw_atom_.find( 'u' );
		} else {
			tr << "WARNING! Asked for vdw_atom_list for " << which_nucleotide << " and it did not exist! " << std::endl;
			utility::vector1< std::string > blank_vector;
			return blank_vector;
		}
	}

	return ( iter->second );
}

Size
rna_residue_type_to_num( chemical::ResidueType const & rt ) {
	if ( rt.name1() == 'a' ) return 1;
	if ( rt.name1() == 'c' ) return 2;
	if ( rt.name1() == 'g' ) return 3;
	if ( rt.name1() == 'u' ) return 4;
	if ( rt.name1() == 't' ) return 4;
	if ( rt.name1() == 'Z' ) return 5; // Mg(2+)
	
	if ( rt.na_analogue() == chemical::na_rad ) return 1;
	if ( rt.na_analogue() == chemical::na_rcy ) return 2;
	if ( rt.na_analogue() == chemical::na_rgu ) return 3;
	if ( rt.na_analogue() == chemical::na_ura ) return 4;
	
	tr << "What is this? " << rt.name() << std::endl;
	utility_exit_with_message( "Asked for rna_residue_name_to_num for unknown residue_name" );
	return 0;
}

Real
RNA_AtomVDW::bump_parameter( Size const atom1, Size const atom2,
	chemical::ResidueType const & rsd1, chemical::ResidueType const & rsd2 ) const
{
	return rna_vdw_parameter_( atom1, atom2,
		rna_residue_type_to_num( rsd1 ),
		rna_residue_type_to_num( rsd2 ) );
}


} //rna
} //scoring
} //core

