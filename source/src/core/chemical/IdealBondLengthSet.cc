// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/IdealBondLengthSet.cc
/// @brief
/// @author Gordon Lemmon

// Unit headers
#include <core/chemical/IdealBondLengthSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

// Project headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>

#include <utility/exit.hh>

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer tr( "core.chemical" );

IdealBondLengthSet::IdealBondLengthSet():
	bond_lengths_()
{}

IdealBondLengthSet::~IdealBondLengthSet() = default;

bool IdealBondLengthSet::contains_bond_length(std::string const & atom_type_name1, std::string const & atom_type_name2) const{
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCOP atom_type_set= cm->atom_type_set( core::chemical::FA_STANDARD );

	AtomTypeIndex atom_type_index1 = atom_type_set->atom_type_index(atom_type_name1);
	AtomTypeIndex atom_type_index2 = atom_type_set->atom_type_index(atom_type_name2);

	return contains_bond_length(atom_type_index1, atom_type_index2);
}

bool IdealBondLengthSet::contains_bond_length(AtomTypeIndex atom_type_index1, AtomTypeIndex atom_type_index2) const{
	if ( atom_type_index1 > atom_type_index2 ) std::swap(atom_type_index1, atom_type_index2);
	std::pair<AtomTypeIndex,AtomTypeIndex> index_pair(atom_type_index1,atom_type_index2);

	if ( bond_lengths_.find(index_pair) == bond_lengths_.end() ) {
		return false;
	}

	return true;
}

BondLength
IdealBondLengthSet::get_bond_length( std::string const & atom_type_name1, std::string const & atom_type_name2) const{
	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCOP atom_type_set= cm->atom_type_set( core::chemical::FA_STANDARD );

	Size atom_type_index1 = atom_type_set->atom_type_index(atom_type_name1);
	Size atom_type_index2 = atom_type_set->atom_type_index(atom_type_name2);

	if ( atom_type_index1 > atom_type_index2 ) std::swap(atom_type_index1, atom_type_index2);

	if ( ! contains_bond_length(atom_type_index1, atom_type_index2) ) {
		utility_exit_with_message( "ideal bond_length not defined between these atom types: "+ atom_type_name1+" "+atom_type_name2);
	}
	return get_bond_length(atom_type_index1, atom_type_index2);
}

BondLength
IdealBondLengthSet::get_bond_length(AtomTypeIndex const atom_type_index1, AtomTypeIndex const atom_type_index2) const{
	debug_assert(contains_bond_length(atom_type_index1, atom_type_index2));

	std::pair<AtomTypeIndex,AtomTypeIndex> index_pair(atom_type_index1,atom_type_index2);

	return bond_lengths_.find(index_pair)->second;
}

void
IdealBondLengthSet::add_bond_length(
	std::string const & atom_type_name1,
	std::string const & atom_type_name2,
	BondLength const length
){

	core::chemical::ChemicalManager *cm= core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCOP atom_type_set= cm->atom_type_set( core::chemical::FA_STANDARD );

	Size atom_type_index1 = atom_type_set->atom_type_index(atom_type_name1);
	Size atom_type_index2 = atom_type_set->atom_type_index(atom_type_name2);

	if ( atom_type_index1 > atom_type_index2 ) std::swap(atom_type_index1, atom_type_index2);

	if ( contains_bond_length(atom_type_index1, atom_type_index2) ) {
		utility_exit_with_message("this pair is already in the table... "+ atom_type_name1+" "+atom_type_name2);
	} else {
		add_bond_length(atom_type_index1, atom_type_index2, length);
	}
}

void
IdealBondLengthSet::add_bond_length(
	AtomTypeIndex const atom_type_index1,
	AtomTypeIndex const atom_type_index2,
	BondLength const length
){
	debug_assert(atom_type_index1 <= atom_type_index2);
	debug_assert(!contains_bond_length(atom_type_index1, atom_type_index2) );

	std::pair<AtomTypeIndex,AtomTypeIndex> index_pair(atom_type_index1,atom_type_index2);
	bond_lengths_[index_pair]= length;
}


void
IdealBondLengthSet::read_file( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() );

	if ( !data.good() ) utility_exit_with_message( "Unable to open element file: "+filename );

	// now parse the rest of the file

	using namespace basic;

	std::string line;
	// parse the header line
	getline( data, line ); // throw out the header line (currently it is just for file readability)
	while ( getline( data,line ) ) {
		utility::trim(line, " \t\n"); // remove leading and trailing spaces
		if ( line.empty() > 0 ) continue; //skip blank lines
		if ( line.find("#",0) == 0 ) continue; // skip comment lines

		std::istringstream l( line );
		std::string name1, name2;
		Real length;

		l >> name1 >> name2 >> length;

		if ( l.fail() ) {
			utility_exit_with_message("bad line: "+line);
		}

		add_bond_length( name1, name2, length);
	}
}


/// @details This function iterates over each bond_type pair in the bond_length_ map and
/// prints both keys. It is only used for debugging.
void
IdealBondLengthSet::print_all_bond_lengths()
{
	for (
			std::map< std::pair<AtomTypeIndex, AtomTypeIndex>, BondLength >::const_iterator i = bond_lengths_.begin(),
			e = bond_lengths_.end();
			i != e;
			++i
			) {
		std::cout << i->first.first << " " << i->first.second << " " << i->second << std::endl;
	}
}

} // chemical
} // core
