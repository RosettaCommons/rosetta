// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/SmartSewingResidue.hh
/// @brief a minimal container for SEWING residues
/// @author frankdt (frankdt@email.unc.edu)


#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <core/conformation/Atom.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
namespace protocols {
namespace sewing {
namespace data_storage {

///@brief a minimal container for SEWING residues


LigandResidue::LigandResidue():
	SmartSewingResidue()
{
	current_contacts_.clear();
	owner_segment_ = nullptr;
	alignment_atoms_.clear();
}
LigandResidue::LigandResidue(LigandResidue const & src):
	SmartSewingResidue( src )
{
	//Make a deep copy
	current_contacts_.clear();
	//current_contacts_ = src.get_current_contacts();
	for ( LigandContactOP contact: src.get_current_contacts() ) {
		current_contacts_.push_back( LigandContactOP( new LigandContact( contact->segment_id, contact->residue_number, contact->residue_atom, contact->ligand_atom ) ) );
	}
	owner_segment_ = src.get_owner_segment();
	ligand_id_ = src.get_ligand_id();
	//TODO
	ideal_contacts_ = src.get_ideal_contacts(); //This one might be OK since they aren't pointers
	alignment_atoms_ = src.get_alignment_atoms();
}

LigandResidue::~LigandResidue()
{}

LigandResidueOP
LigandResidue::clone() const{
	return LigandResidueOP( new LigandResidue( *this ) );
}

//Getters
utility::vector1< LigandContactOP > const &
LigandResidue::get_current_contacts() const{
	return current_contacts_;
}

utility::vector1< LigandContactOP > &
LigandResidue::get_nonconst_current_contacts(){
	return current_contacts_;
}
/* utility::vector1< LigandContact >
LigandResidue::get_preferred_contacts(){
return preferred_contacts_;
}
*/

data_storage::LigandSegmentOP
LigandResidue::get_owner_segment() const {
	return owner_segment_;
}

data_storage::LigandSegmentOP
LigandResidue::get_nonconst_owner_segment() {
	return owner_segment_;
}

//Setters

void
LigandResidue::set_owner_segment( data_storage::LigandSegmentOP owner ){
	owner_segment_ = owner;
}

void
LigandResidue::add_contact( LigandContactOP new_contact ) {
	current_contacts_.push_back( new_contact );
}

void
LigandResidue::set_contacts( utility::vector1< LigandContactOP > contacts ){
	current_contacts_ = contacts;
}

core::Size
LigandResidue::get_ligand_id() const{
	return ligand_id_;
}
void
LigandResidue::set_ligand_id( core::Size id ){
	ligand_id_ = id;
}



void
LigandResidue::set_alignment_atoms( utility::vector1< core::Size > atoms ){
	alignment_atoms_ = atoms;
}

utility::vector1< core::Size >
LigandResidue::get_alignment_atoms() const{
	return alignment_atoms_;
}


void
LigandResidue::add_ideal_contact( hashing::IdealContact contact ){
	core::Size contact_atom = contact.ligand_atom;
	ideal_contacts_[ contact_atom ] = contact;
}

std::map< core::Size, hashing::IdealContact > const &
LigandResidue::get_ideal_contacts() const{
	return ideal_contacts_;
}
/*
void
LigandResidue::add_preferred_contact( LigandContact contact ){
preferred_contacts_.push_back( contact );
}
*/
} //protocols
} //sewing
} //data_storage





