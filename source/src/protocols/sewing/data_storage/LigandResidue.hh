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


#ifndef INCLUDED_protocols_sewing_data_storage_LigandResidue_hh
#define INCLUDED_protocols_sewing_data_storage_LigandResidue_hh

#include <protocols/sewing/data_storage/LigandResidue.fwd.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/LigandSegment.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/conformation/Residue.fwd.hh>
namespace protocols {
namespace sewing {
namespace data_storage {

///@brief a minimal container for SEWING residues
class LigandResidue : public SmartSewingResidue {

public:

	LigandResidue();
	LigandResidue(LigandResidue const & src);

	virtual ~LigandResidue();

	LigandResidueOP
	clone() const;

	//Getters
	utility::vector1< LigandContactOP >  const &
	get_current_contacts() const;

	utility::vector1< LigandContactOP > &
	get_nonconst_current_contacts();
	data_storage::LigandSegmentOP
	get_owner_segment() const;

	data_storage::LigandSegmentOP
	get_nonconst_owner_segment();
	/*
	utility::vector1< LigandContact >
	get_preferred_contacts();
	*/
	//Setters
	void
	add_contact( LigandContactOP contact );

	void
	set_contacts( utility::vector1< LigandContactOP > contacts );

	void
	set_owner_segment( data_storage::LigandSegmentOP owner );
	//void
	//add_preferred_contact( LigandContact );

	core::Size
	get_ligand_id() const;
	void
	set_ligand_id( core::Size );

	utility::vector1< core::Size >
	get_alignment_atoms() const;

	void
	set_alignment_atoms( utility::vector1< core::Size > );

	void
	add_ideal_contact( hashing::IdealContact contact );

	std::map< core::Size, hashing::IdealContact > const &
	get_ideal_contacts() const;

	bool
	get_partner_ligand() const{ return partner_ligand_; }

	void
	set_partner_ligand( bool input ){ partner_ligand_ = input; }


private:
	//utility::vector1< data_storage::SmartSewingResidueCOP > current_contacts_;
	utility::vector1< LigandContactOP > current_contacts_; //Pairs are (segment ID, residue number )
	data_storage::LigandSegmentOP owner_segment_;
	core::Size ligand_id_;
	utility::vector1< core::Size > alignment_atoms_;
	std::map< core::Size, hashing::IdealContact > ideal_contacts_;
	bool partner_ligand_ = false;
};

struct LigandContact
{
	core::Size ligand_atom; //By default anchor to first atom in ligand
	core::Size residue_atom; //By default anchor to CA
	core::Size residue_number; //You have to supply a residue number
	core::Size segment_id; //You have to supply a segment ID
	LigandContact( core::Size segid, core::Size resnum, core::Size resat, core::Size ligat ){
		ligand_atom = ligat;
		segment_id = segid;
		residue_atom = resat;
		residue_number = resnum;
	}
	LigandContact( core::Size segid, core::Size resnum ){
		segment_id = segid;
		residue_number = resnum;
	}
};

//This data structure is used for describing ligand contacts before they are made!
struct ContactDescription
{
	bool partner_contact=false;
	std::string contact_resnum_string="";
	std::string ligand_atom_name="";
	std::string contact_atom_name="";
	core::Size contact_resnum=0;
	core::Size contact_atom_num=0;
	core::Size ligand_atom_num=0;
};

struct LigandDescription
{
	core::Size ligand_id = 0;
	bool partner_ligand=false;
	std::string ligand_resnum_string="";
	core::select::residue_selector::ResidueSelectorCOP ligand_selector=nullptr;
	core::Size ligand_resnum=0;
	bool auto_detect_contacts=true;
	utility::vector1< ContactDescription > ligand_contacts;
	utility::vector1< std::string > alignment_atoms_str;
	utility::vector1< core::Size > alignment_atoms_num;
	std::map< std::string, hashing::IdealContact > ideal_contacts_str;
	std::map< core::Size, hashing::IdealContact > ideal_contacts_num;
	std::string pdb_conformers_string;
	utility::vector1< core::conformation::ResidueCOP > pdb_conformers;
	utility::vector1< std::string > ligand_coord_files;
	utility::pointer::shared_ptr< std::map< char, hashing::LigandHashMap > > ligand_coords;
	core::Real geometry_score_threshold;
	LigandDescription(){}
	LigandDescription(LigandDescription const & other ):
		ligand_id( other.ligand_id ),
		partner_ligand( other.partner_ligand ),
		ligand_resnum_string( other.ligand_resnum_string ),
		ligand_selector( other.ligand_selector ),
		ligand_resnum( other.ligand_resnum ),
		auto_detect_contacts( other.auto_detect_contacts ),
		ligand_contacts( other.ligand_contacts ),
		alignment_atoms_str( other.alignment_atoms_str ),
		alignment_atoms_num( other.alignment_atoms_num ),
		ideal_contacts_str( other.ideal_contacts_str ),
		ideal_contacts_num( other.ideal_contacts_num ),
		pdb_conformers_string( other.pdb_conformers_string ),
		pdb_conformers( other.pdb_conformers ),
		ligand_coord_files( other.ligand_coord_files ),
		ligand_coords( other.ligand_coords ),
		geometry_score_threshold( other.geometry_score_threshold )
	{
	}
};



} //protocols
} //sewing
} //data_storage



#endif //INCLUDED_protocols_sewing_data_storage_LigandResidue_hh





