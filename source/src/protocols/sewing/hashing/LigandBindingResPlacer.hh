// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LigandBindingResPlacer.hh
/// @brief Using the file format set in zinc_statistic_generator, checks for residues that could coordinate a metal.
/// @author guffysl (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_hashing_LigandBindingResPlacer_hh
#define INCLUDED_protocols_sewing_hashing_LigandBindingResPlacer_hh

#include <protocols/sewing/hashing/LigandBindingResPlacer.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>
#include <protocols/sewing/data_storage/LigandResidue.fwd.hh>
#include <protocols/sewing/data_storage/LigandSegment.fwd.hh>
//#include <protocols/sewing/data_storage/SmartSewingResidue.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/types.hh>

#include <numeric/xyzVector.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>




namespace protocols {
namespace sewing {
namespace hashing {


struct LigandCoordInfo {
	numeric::xyzVector< core::Real > local_coords;
	std::string coord_res_name;
	core::Size coord_atom_index;
	core::Size ligand_atom_index;
	utility::fixedsizearray1< core::Real, 4 > chis;


	friend
	bool
	operator< (LigandCoordInfo const & a, LigandCoordInfo const & b ){

		if ( a.ligand_atom_index < b.ligand_atom_index ) return true;
		if ( b.ligand_atom_index < a.ligand_atom_index ) return false;

		if ( a.coord_atom_index < b.coord_atom_index ) return true;
		if ( b.coord_atom_index < a.coord_atom_index ) return false;
		//The following is for coordination between the same atom pair

		if ( a.local_coords.x() < b.local_coords.x() ) return true;
		if ( b.local_coords.x() < a.local_coords.x() ) return false;

		if ( a.local_coords.y() < b.local_coords.y() ) return true;
		if ( b.local_coords.y() < a.local_coords.y() ) return false;

		if ( a.local_coords.z() < b.local_coords.z() ) return true;
		if ( b.local_coords.z() < a.local_coords.z() ) return false;

		if ( a.chis[ 1 ] < b.chis[ 1 ] ) return true;
		if ( b.chis[ 1 ] < a.chis[ 1 ] ) return false;

		if ( a.chis[ 2 ] < b.chis[ 2 ] ) return true;
		if ( b.chis[ 2 ] < a.chis[ 2 ] ) return false;

		if ( a.chis[ 3 ] < b.chis[ 3 ] ) return true;
		if ( b.chis[ 3 ] < a.chis[ 3 ] ) return false;

		if ( a.chis[ 4 ] < b.chis[ 4 ] ) return true;
		if ( b.chis[ 4 ] < a.chis[ 4 ] ) return false;
		return false;
	}

};



typedef boost::unordered_map< HashKey, std::set< LigandCoordInfo >, coord_hash, coord_equal_to > LigandHashMap;

///@details Given a ligand and an assembly, this class identifies the potential best mutations (with chi angles) to make to add a coordinating residue to the ligand. Another class will test the geometries and make the mutations.
class LigandBindingResPlacer: public utility::pointer::ReferenceCount {

public:

	///@brief Default constructor
	LigandBindingResPlacer();

	///@brief Constructor with arguments
	LigandBindingResPlacer( utility::pointer::shared_ptr< std::map< char, LigandHashMap > > coords, data_storage::LigandResidueOP ligand_residue );

	///@brief Parses the ligand file format
	static void
	parse_ligand_files( utility::vector1< data_storage::LigandDescription > & ligands );
	///@details
	//The form it should be parsed into should fill the following criteria:
	//1) Organized by secondary structure
	//2) Organized by residue to avoid unnecessary checks
	//3) Ligand coordinates hashed for fast comparison using numeric::hash_value( xyzVector ) for first step; stored as vector for second step
	//4) We won't need chi angles until the second step

	///First format: Map of (secstruct, map( restype, vector< LigandCoordInfo> ) )
	///Second format: Map of secstruct to a std::set of zinc HashKeys (generated the same basic way we do in Hasher)
	///@brief Checks all residue stubs in assembly vs metal coordinates to determine which residues are potentially metal-coordinating. Return a set of SmartSewingResidueOP.
	//std::set< data_storage::SmartSewingResidueOP >
	std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> >
	identify_possible_binders( data_storage::SmartAssemblyCOP assembly );
	///@details Iterate over all segments
	///Iterate over all residues in the segment
	///If it's a required residue, skip it
	///Otherwise, convert it to a stub
	///Transform the ligand's coordinates into the stub's local coordinates

	//How do we deal with multi-atom ligands?
	//We'll need to check all the atoms against the LigandHashMap


	///Convert the ligand atoms' coordinates into HashKeys
	///Look for the HashKeys in the appropriate HashMap for this DSSP, and compare the atom name to see if it's a real match
	///If it's there, add this residue to the map with the corresponding entries from LigandCoordInfo
	//The core::Size is the residue's position within the Assembly--add a function to SmartAssembly that get the nth residue in the assembly



	core::conformation::ResidueOP
	make_aligned_residue( data_storage::SmartSewingResidueOP, std::string );

	///@brief Finds the best residue/rotamer to coordinate metal (Main function for class)
	std::pair< data_storage::SmartSegmentOP, bool >
	choose_best_metal_coordinator( data_storage::SmartAssemblyOP assembly );
	///@details First identify possible metal binding residues
	///Iterate over all residues in the set
	///Find best coordinating rotamer for each residue per restype


	///@brief Generates a HashKey from an xyzVector of core::Real
	static
	HashKey
	generate_ligand_key( numeric::xyzVector< core::Real > );
	///@brief Finds the best residue/rotamer to coordinate metal for a given residue and returns for each residue type a pair of (LigandCoordInfo, distance) where distance is the distance of the metal from its ideal position
	std::pair< LigandCoordInfo, core::Real >
	best_rotamer_for_residue( std::pair< core::Size, core::Size > res, std::set< LigandCoordInfo > res_coord, data_storage::SmartAssemblyOP assembly );
	//At this point we've already created ligand_local_coords. We'll need to also know the DSSP of this residue, so maybe we should just take a SmartSegment and residue number?
	///@details For the given residue, does the following:
	///1) Looks the residue up in possible_coordinations_
	///Make a residue of the appropriate type and transform it onto res
	///Iterate over LigandCoorInfo and for each one set the chi angles



	void
	replace_chimaera_ligand_parent( data_storage::SmartAssemblyOP assembly, data_storage::SmartSegmentOP chimaeric_segment, data_storage::LigandSegmentOP replacement_ligseg,  bool change_n_terminal_parent, bool first_basis_is_n_terminal, core::Size chimaeric_contact_resnum, utility::vector1< std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > > > res_scores );


	void
	revert_last_added_contact( data_storage::SmartAssemblyOP assembly );

	//Getters
	data_storage::LigandResidueOP
	get_ligand() const;

	std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> > const &
	get_possible_coordinations() const;

	utility::pointer::shared_ptr< std::map< char, LigandHashMap > const >
	get_ligand_local_coords() const;

	core::Size
	contacts_remaining() const;

	core::Real
	get_geometry_score_weight() const;

	core::Real
	get_best_geometry_score() const;

	//Setters
	void
	set_ligand( data_storage::LigandResidueOP lig );

	void
	set_geometry_score_weight( core::Real );

	void
	set_best_geometry_score( core::Real );

	void
	set_ligand_local_coords( utility::pointer::shared_ptr< std::map< char, LigandHashMap > > ligand_local_coords );

	//Static functions for xml
	static void
	parse_ideal_contacts( utility::tag::TagCOP tag, utility::vector1< data_storage::LigandDescription > & ligands );

	static void
	define_ideal_contacts_subelement( utility::tag::XMLSchemaSimpleSubelementList & subs, utility::tag::XMLSchemaDefinition & xsd );

	static std::string
	ligand_subtag_ct_namer( std::string tag );
private:

	//Information I will need
	//LigandResidueOP
	data_storage::LigandResidueOP ligand_;
	//Skip any required residues (which include residues already bound to a ligand)
	//Data formats created in parsing method
	//Second format: Map of secstruct to a std::set of zinc HashKeys (generated the same basic way we do in Hasher)
	utility::pointer::shared_ptr< std::map< char, LigandHashMap > > ligand_local_coords_;
	//Information on ideal ligand coordination
	/*
	std::map< core::Size, IdealContact > ideal_contacts_;
	*/
	std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> > possible_coordinations_;
	core::Size contacts_to_make_;

	data_storage::LigandContactOP last_contact_added_;
	core::Real geometry_score_weight_;
	core::Real best_geometry_score_ = 0;

	core::chemical::AtomTypeSetCAP atom_types_;


};




} //namespace hashing
} //namespace sewing
} //namespace protocols

#endif //include guard
