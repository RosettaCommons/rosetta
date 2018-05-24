// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/SmartSewingResidue.cc
/// @brief a minimal container for SEWING residues
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sewing.data_storage.SmartSewingResidue" );


namespace protocols {
namespace sewing {
namespace data_storage {

SmartSewingResidue::SmartSewingResidue():
	utility::pointer::ReferenceCount()
{
	chis_.clear();
	amino_acid_type_="NNN";
	full_type_name_="NNN";
}

SmartSewingResidue::~SmartSewingResidue(){}

SmartSewingResidue::SmartSewingResidue( SmartSewingResidue const & other ) {
	//Deep copy the atom vector
	atom_vector_.clear();
	for ( core::conformation::Atom atom: other.get_const_atom_vector() ) {
		//We're going to try making a whole new xyzVector
		numeric::xyzVector< core::Real > new_coords;
		new_coords.x( atom.xyz().x() );
		new_coords.y( atom.xyz().y() );
		new_coords.z( atom.xyz().z() );
		atom_vector_.push_back( core::conformation::Atom( new_coords, atom.type(), atom.mm_type() ) );
	}
	chis_ = other.get_chis();
	amino_acid_type_ = other.get_amino_acid_type();
	if ( other.get_full_type_name() == "NNN" ) { //For use with old segment files
		full_type_name_ = amino_acid_type_;
	} else {
		full_type_name_ = other.get_full_type_name();
	}

}



SmartSewingResidueOP
SmartSewingResidue::clone() const {
	return SmartSewingResidueOP( new SmartSewingResidue( *this ) );
}

void
SmartSewingResidue::set_atom_vector(utility::vector1<core::conformation::Atom> new_atom_vector) {
	atom_vector_ = new_atom_vector; // there should be an "apply transform" function instead of this
}

utility::vector1<core::conformation::Atom> &
SmartSewingResidue::get_atom_vector(){
	return atom_vector_;
}
//const version
utility::vector1<core::conformation::Atom> const &
SmartSewingResidue::get_const_atom_vector() const{
	return atom_vector_;
}

core::conformation::Atom &
SmartSewingResidue::get_atom(core::Size atom_number){
	return atom_vector_[atom_number];
}
void
SmartSewingResidue::set_chis(utility::vector1<core::Real> new_chis){
	chis_ = new_chis;
}

utility::vector1<core::Real>
SmartSewingResidue::get_chis() const{
	return chis_;
}

void
SmartSewingResidue::set_amino_acid_type(std::string new_amino_acid_type){
	amino_acid_type_ = new_amino_acid_type;
}

std::string
SmartSewingResidue::get_amino_acid_type() const{
	return amino_acid_type_;
}

void
SmartSewingResidue::become(SmartSewingResidueCOP residue_to_become){
	if ( residue_to_become == nullptr ) {
		chis_.clear();
		amino_acid_type_ = "";
		full_type_name_ = "";
		atom_vector_.clear();
	} else {
		//TR << this->get_amino_acid_type() << std::endl;
		//TR << this->get_chis() << std::endl;
		//TR << "Becoming Residue" << std::endl;
		//TR << residue_to_become->get_chis() << std::endl;
		//TR << residue_to_become->get_amino_acid_type() << std::endl;
		this->set_chis(residue_to_become->get_chis());
		this->set_amino_acid_type(residue_to_become->get_amino_acid_type());
		this->set_full_type_name( residue_to_become->get_full_type_name() );
		this->get_atom_vector().clear();
		//Deep copy the atom vector as well
		for ( core::conformation::Atom old_atom: residue_to_become->get_const_atom_vector() ) {
			core::conformation::Atom atom_copy( old_atom );
			this->get_atom_vector().push_back( atom_copy );
		}
		this->set_atom_vector(residue_to_become->get_const_atom_vector());
	}
}

void
SmartSewingResidue::set_type( PolymericType type ){
	type_ = type;
}

PolymericType
SmartSewingResidue::get_type(){
	return type_;
}

std::string
SmartSewingResidue::get_full_type_name() const{
	return full_type_name_;
}

void
SmartSewingResidue::set_full_type_name( std::string const & type ){
	full_type_name_ = type;
}

} //protocols
} //sewing
} //data_storage






