// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/LigandSegment.hh
/// @brief a neighbor-aware SewSegment version
/// @author frankdt (frankdt@email.unc.edu)


#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>

namespace protocols {
namespace sewing {
namespace data_storage {
///@brief a segment that contains a single ligand residue with given attachment point(s)

LigandSegment::LigandSegment():// "everything is empty" default constructor. This gets used for building segments from files
	SmartSegment()
{
	ligand_residues_.clear();
	ligand_contacts_.clear();
	owned_ligand_residues_.clear();
}
LigandSegment::LigandSegment( bool is_vital ):
	SmartSegment( is_vital )
{
	ligand_residues_.clear();
	ligand_contacts_.clear();
	owned_ligand_residues_.clear();
}
LigandSegment::LigandSegment(bool is_vital, core::Size max_segment_length): // use this constructor to actually make segments in an assembly, since it preallocates everything
	SmartSegment( is_vital, max_segment_length )
{
	ligand_residues_.clear();
	ligand_contacts_.clear();
	owned_ligand_residues_.clear();
}

LigandSegment::LigandSegment(LigandSegment const & src):
	SmartSegment( src )
{
	std::map< core::Size, LigandResidueOP > owned_ligs = src.get_const_owned_ligand_residues();
	for ( std::pair< core::Size, LigandResidueOP > owned_ligand: owned_ligs ) {
		LigandResidueOP new_owned_ligand = owned_ligand.second->clone();
		attach_ligand( new_owned_ligand, true );
	}
	//In this case, we don't have pointers to unowned ligands
	//Since we're cloning, we can assume these ligands exist and just add their ids
	for ( core::Size lig: src.get_ligand_residues() ) {
		ligand_residues_.insert( lig );
	}
	for ( core::Size contact: src.get_ligand_contact_indices() ) {
		add_ligand_contact( contact );
	}
}


///@details This constructor is for use with LigandBindingResPlacer
///It knows:
///*pointer to ligand (and therefore the ligand ID)
///*ligand contact residue(because it just added it)
///*There won't be any owned ligands
///We only need to allow for one ligand since this will only be used when the first contact is added
LigandSegment::LigandSegment(SmartSegment const & src, bool is_vital, std::set< core::Size > ligand_binding_residues, core::Size ligand_to_attach):
	SmartSegment( src )
{
	this->set_is_vital( is_vital );
	for ( core::Size res: ligand_binding_residues ) {
		add_ligand_contact( res );
	}
	ligand_residues_.insert( ligand_to_attach );
}


LigandSegment::~LigandSegment(){
	//this->~SmartSegment();
	//Maybe set owned ligands' owner segments to nullptr, but maybe not (if it's possible more than one seg thinks it owns a ligand)
}

LigandSegmentOP
LigandSegment::clone() const{
	return LigandSegmentOP( new LigandSegment( *this ) );
}

std::set< core::Size >
LigandSegment::get_ligand_residues() const{
	return ligand_residues_;
}

std::set< core::Size > &
LigandSegment::get_nonconst_ligand_residues(){
	return ligand_residues_;
}

std::set< core::Size >
LigandSegment::get_ligand_contact_indices() const{
	return ligand_contacts_;
}

std::set< core::Size > &
LigandSegment::get_nonconst_ligand_contact_indices(){
	return ligand_contacts_;
}
void
LigandSegment::add_ligand_contact( core::Size contact_index ){
	ligand_contacts_.insert( contact_index );
	this->add_vital_residue( contact_index );
}

void
LigandSegment::attach_ligand( LigandResidueOP ligand_res, bool owner ){
	//Add the ligand residue
	if ( owner ) {
		owned_ligand_residues_[ ligand_res->get_ligand_id() ] = ligand_res;
	}
	//Unowned ligands should only be replaced if they aren't already in there
	ligand_residues_.insert( ligand_res->get_ligand_id() );
}


void
LigandSegment::attach_unowned_ligand( LigandResidueCOP ligand_res ){ //We still require a pointer to the ligand in order to attach it--prevents attaching nonexistent ligands
	ligand_residues_.insert( ligand_res->get_ligand_id() );
}


std::map< core::Size, LigandResidueOP > &
LigandSegment::get_owned_ligand_residues(){
	return owned_ligand_residues_;
}

std::map< core::Size, LigandResidueOP > const
LigandSegment::get_const_owned_ligand_residues() const{
	return owned_ligand_residues_;
}

} //protocols
} //sewing
} //data_storage







