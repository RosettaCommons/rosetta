// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/AtomVDW.hh
/// @brief
/// @author Phil Bradley
/// @author Rhiju Das


// Unit Headers
#include <core/scoring/rna/RNA_AtomVDW.hh>

// Package headers

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace rna {

Size const
rna_residue_name_to_num( char const c )
{
	if ( c == 'a') return 1;
	if ( c == 'c') return 2;
	if ( c == 'g') return 3;
	if ( c == 'u') return 4;
	return 0;
}

//Why doesn't this helper function already exist in vector class?
Size
get_position_in_vector( utility::vector1< std::string> & vec, std::string const element )
{
	Size count( 1 );
	for ( utility::vector1<std::string>::iterator iter = vec.begin(); iter < vec.end(); iter++ )	{
		if (*iter == element) return count;
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
	basic::database::open( stream, "rna_atom_vdw.txt" );

	rna_vdw_parameter_.dimension( 10, 10, 4, 4 );

	if ( !stream.good() ) utility_exit_with_message( "Unable to open rna_atom_vdw.txt!" );

	// read the entire file and figure out what atom_types are present and in what order
	utility::vector1< std::string > lines;

	std::string line;
	std::string atom_name1, atom_name2;
	char which_residue1, which_residue2;
	Real input_bump_parameter;
	while ( getline( stream, line ) ) {
		lines.push_back(line);
		std::istringstream l(line);
		l >> which_residue1 >> which_residue2 >> atom_name1 >> atom_name2 >> input_bump_parameter;

		Size const pos1 = get_position_in_vector( rna_vdw_atom_[which_residue1], atom_name1 );
		Size const pos2 = get_position_in_vector( rna_vdw_atom_[which_residue2], atom_name2 );

		rna_vdw_parameter_( pos1,
												pos2,
												rna_residue_name_to_num( which_residue1 ),
												rna_residue_name_to_num( which_residue2 ) ) = input_bump_parameter;

	}

	num_rna_vdw_atoms_check_ = rna_vdw_atom_['g'].size();

	initialize_atom_numbers(); // For fast look up of atom indices

}

///////////////////////////////////////////////
void
RNA_AtomVDW::initialize_atom_numbers()
{

	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	initialize_atom_numbers( rsd_set, na_rad );
	initialize_atom_numbers( rsd_set, na_rgu );
	initialize_atom_numbers( rsd_set, na_rcy );
	initialize_atom_numbers( rsd_set, na_ura );

}

///////////////////////////////////////////////
void
RNA_AtomVDW::initialize_atom_numbers( chemical::ResidueTypeSetCAP & rsd_set, chemical::AA const & aa )
{
	using namespace core::chemical;

	ResidueTypeCAPs const & rsd_types( rsd_set->aa_map( aa ) );

	if ( rsd_types.size() < 1) {
		std::cout << "PROBLEM FINDING RESIDUE: " << name_from_aa( aa ) << std::endl;
		return;
	}

	ResidueTypeCAP const & rsd_type( rsd_types[1] );
	char const which_nucleotide( rsd_type->name1() ) ;

	utility::vector1 < Size > atom_numbers_temp;

	for (Size i = 1; i <= num_rna_vdw_atoms_check_; i++ ) {
		atom_numbers_temp.push_back( rsd_type->atom_index( rna_vdw_atom_[ which_nucleotide ][i] ) );
	}

	atom_numbers_[ which_nucleotide ] =  atom_numbers_temp;

}


// std::string
// RNA_AtomVDW::check_atom( Size const atom_index, char const which_nucleotide ) const
// {
// 	// Following doesn't work, unfortunately -- not actually a const operation for maps.
// 	//	return rna_vdw_atom_[ which_nucleotide ][ atom_index ];

// 	AtomList::const_iterator iter = rna_vdw_atom_.find( which_nucleotide );
//   return (iter->second)[ atom_index ];
// }

//////////////////////////////////////////////
utility::vector1 < std::string > const
RNA_AtomVDW::vdw_atom_list( char const which_nucleotide ) const
{
	AtomList::const_iterator iter = rna_vdw_atom_.find( which_nucleotide );
	return (iter->second);
}

Real
RNA_AtomVDW::bump_parameter( Size const atom1, Size const atom2,
														 char const which_residue1, char const which_residue2 ) const
{
	return rna_vdw_parameter_( atom1, atom2,
														 rna_residue_name_to_num( which_residue1 ),
														 rna_residue_name_to_num( which_residue2 ) );
}

//////////////////////////////////////////////
utility::vector1< Size > const &
RNA_AtomVDW::atom_numbers_for_vdw_calculation( char const which_nucleotide ) const
{
	AtomNumberList::const_iterator iter = atom_numbers_.find( which_nucleotide );
	return (iter->second);
}


} // namespace rna
} // namespace scoring
} // namespace core

