// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Atom.cc
/// @brief  Energy graph class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/chemical/Atom.hh>

// Boost Headers
#include <core/graph/unordered_object_pool.hpp>

#include <iostream>

#include <utility/vector1.hh>
#include <boost/pool/pool.hpp>


namespace core {
namespace chemical {

using namespace graph;

///////// Atom Class /////////////

Atom::Atom( Graph * owner, Size index ) :
	parent( owner, index )//, moved_( false )
{}

Atom::Atom(
		Graph * owner,
		Size index,
		std::string const & name_in,
	//	std::string const type_name,
		std::string const mm_name,
		Size const atom_type_index,
		Size const mm_atom_type_index,
		Real const charge,
		Vector const ideal_xyz,
		AtomICoor const icoor

):
	name_( name_in ),
	//type_name_(type_name),
	mm_name_(mm_name),
	atom_type_index_(atom_type_index),
	mm_atom_type_index_(mm_atom_type_index),
	charge_(charge),
	ideal_xyz_(ideal_xyz),
	icoor_(icoor)
{}

Atom::~Atom() {}


void Atom::print() const
{
	std::cout << "Atom::print() deferring to parent::print()" << std::endl;
	parent::print();
}

/// @brief copy mmember data from source node
///
/// invoked by copy ctor and operator= methods from Graph base class
void Atom::copy_from( parent const * source )
{
	//Atom const * en_source = utility::down_cast< Atom const * > ( source ); //nothing to copy, still -- want to assert the dynamic cast
	Atom const * atom = static_cast< Atom const * > ( source );
	name_ = atom->name_ ;
	//type_name_(atom.type_name),
	mm_name_ = atom->mm_name_;
	atom_type_index_ = atom->atom_type_index_;
	mm_atom_type_index_ =atom->mm_atom_type_index_;
	charge_ = atom->charge_;
	ideal_xyz_ = atom->ideal_xyz_;
	icoor_ = atom->icoor_;
}

Size Atom::count_static_memory() const
{
	return sizeof ( Atom );
}

Size Atom::count_dynamic_memory() const
{
	Size tot = 0;
	tot += parent::count_dynamic_memory(); //recurse to parent
	return tot;
}

//bool Atom::moved() const { return moved_; }
//void Atom::moved( bool setting ) { moved_ = setting; }


} //namespace chemical
} //namespace core
