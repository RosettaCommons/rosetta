// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueTypeSet_hh
#define INCLUDED_core_chemical_ResidueTypeSet_hh


// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>


// Package headers
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
// Commented by inclean daemon #include <core/chemical/AtomTypeSet.hh>
// Commented by inclean daemon #include <core/chemical/MMAtomTypeSet.hh>

// Project headers
//XRW_B_T1
// Commented by inclean daemon #include <core/coarse/TranslatorSet.fwd.hh>
//XRW_E_T1
//#include <core/pose/Pose.fwd.hh>

// Utility Headers
// Commented by inclean daemon #include <utility/pointer/ReferenceCount.hh>
// Commented by inclean daemon #include <utility/vector1_bool.hh>


// STL headers
#include <list>

// AUTO-REMOVED #include <core/chemical/Adduct.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>

//Auto Headers
#include <core/chemical/Adduct.fwd.hh>
#ifdef PYROSETTA
	#include <core/chemical/ResidueType.hh>
#endif



namespace core {
namespace chemical {

///  A collection of ResidueType defined
/**
	 One thing that is not nailed down is whether a single ResidueSet can have ResidueType's with
	 different AtomTypeSets. I've left open this possibility currently although there isnt any code
	 for this yet (PB-07/07)
**/


class ResidueTypeSet : public utility::pointer::ReferenceCount
{
public:
	typedef std::list< AA >::const_iterator AAsIter;
	typedef	std::map< std::string, ResidueTypeCAP >::const_iterator const_residue_iterator;

public:


	/// @brief default c-tor
	ResidueTypeSet();

	/// @brief constructor from directory
	ResidueTypeSet(
		std::string const & name,
		std::string const & directory,
		std::vector< std::string > const & extra_res_param_files = std::vector< std::string >(),
		std::vector< std::string > const & extra_patch_files = std::vector< std::string >()
	);

	virtual ~ResidueTypeSet();

	/// @brief name of the residue type set
	std::string const &
	name() const {
		return name_;
	}

	/// @brief read a list of residue types
	void
	read_list_of_residues(
		std::string const & list_filename,
		AtomTypeSetCAP atom_types,
		ElementSetCAP elements,
		MMAtomTypeSetCAP mm_atom_types,
		orbitals::OrbitalTypeSetCAP orbital_types//,
//		CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd atomtypes until I have had a chance to fully implement them
	);

	void
	read_files(
		utility::vector1< std::string > const & filenames,
		AtomTypeSetCAP atom_types,
		ElementSetCAP elements,
		MMAtomTypeSetCAP mm_atom_types,
		orbitals::OrbitalTypeSetCAP orbital_types//,
//		CSDAtomTypeSetCAP csd_atom_types kwk commenting out the csd atomtypes until I have had a chance to implement them
	);

// Old code - doesn't respect various command line options. See constructor for current patch-loading functionality
//	/// @brief apply patches to base ResidueType to generate variant ResidueTyes
//	void
//	apply_patches( std::string const & filename );

	void
	apply_patches(
		utility::vector1< std::string > const & filenames
	);

	/// @brief apply patches to base ResidueType to generate variant ResidueTyes
	void
	place_adducts();

	/// @brief adds a new residue type to the set
	void
	add_residue_type( ResidueTypeOP new_type );

	/// @brief Create correct combinations of adducts for a residue type
	void create_adduct_combinations(
		ResidueType const & rsd,
		std::map< std::string, int > ref_map,
		std::map< std::string, int > count_map,
		utility::vector1< bool > add_mask,
		utility::vector1< Adduct >::const_iterator work_iter
	);

	/// @brief query ResidueTypes by their 3-letter name
	///
	/// @details 3-letter name is not unique to each ResidueType
	/// for example, 3-letter name "HIS" matches both his tautomers,
	/// HIS and HIS_D. Return an empty list if no match is found.
	ResidueTypeCAPs const &
	name3_map( std::string const & name ) const
	{
		assert( name.size() == 3 );
		if ( name3_map_.find( name ) == name3_map_.end() ) {
			return empty_residue_list_;
		}
		return name3_map_.find( name )->second;
	}

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	ResidueType const &
	name_map( std::string const & name ) const
	{
		if ( name_map_.find( name ) == name_map_.end() ) {
			utility_exit_with_message( "unrecognized residue name '"+name+"'" );
		}
		return *( name_map_.find( name )->second );
	}

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	ResidueType &
	nonconst_name_map( std::string const & name )
	{
		if ( nonconst_name_map_.find( name ) == nonconst_name_map_.end() ) {
			utility_exit_with_message( "unrecognized residue name" );
		}
		return *( nonconst_name_map_.find( name )->second );
	}

	/// @brief query if a ResidueType of the unique residue id is present.
	///
	///
	bool
	has_name( std::string const & name ) const
	{
		return ( name_map_.find( name ) != name_map_.end() );
	}

	/// @brief query if a ResidueType of a certain 3letter name is present.
	bool
	has_name3( std::string const & name3 ) const
	{
		return ( name3_map_.find( name3 ) != name3_map_.end() );
	}

	/// @brief query a variant ResidueType by its base ResidueType and VariantType
	///
	/// @note currently, this will not work for variant types defined as alternate
	/// base residues (ie different params files)
	ResidueType const &
	get_residue_type_with_variant_added( ResidueType const & init_rsd, VariantType const & new_type ) const;


	/// @brief return the residuetype we get from variant rsd type after removing the desired variant type
	///
	/// @note currently, this will not work for variant types defined as alternate
	/// base residues (ie different params files)
	ResidueType const &
	get_residue_type_with_variant_removed( ResidueType const & init_rsd, VariantType const & old_type ) const;


	/// @brief query ResidueTypes by their AA enum type
	///
	/// @details similar to name3_map, return all matched residue types
	/// or an empty list.
	ResidueTypeCAPs const &
	aa_map( AA const & aa ) const
	{
		if ( aa_map_.find( aa ) == aa_map_.end() ) {
			return empty_residue_list_;
		}
		return aa_map_.find( aa )->second;
	}

	/// @brief select a set of ResidueTypes give certain criteria
	void
	select_residues(
		ResidueSelector const & selector,
		ResidueTypeCAPs & matches
	) const;

	/// @brief beginning of aas_defined_ list
	std::list< AA >::const_iterator
	aas_defined_begin() const;

	/// @brief end of aas_defined_ list
	std::list< AA >::const_iterator
	aas_defined_end() const;

	const_residue_iterator
	all_residues_begin() const
	{
		return name_map_.begin();
	}

	const_residue_iterator
	all_residues_end() const
	{
		return name_map_.end();
	}

	/// alternate access to all residuetypes as vector
	ResidueTypeCAPs const &
	residue_types() const
	{
		return residue_types_const_;
	}

	/// @brief accessor for database_directory
	std::string const&
	database_directory() const
	{
		return database_directory_;
	}

	//////////////////
	// private methods
private:

	/// @brief clear residue maps
	void
	clear_residue_maps();

	/// @brief update residue maps
	void
	update_residue_maps();

	void
	add_residue_type_to_maps( ResidueType & rsd );



	//////////////////
	// data
private:

	//chemical::AtomTypeSet const & atom_types_;

// 	/// @brief the atom-types
// 	chemical::AtomTypeSetCAP atom_types_;

// 	/// @brief the MMatom-types
// 	chemical::MMAtomTypeSetCAP mm_atom_types_;

	/// What does the ChemicalManager call this ResidueTypeSet?
	std::string name_;

	/// @brief the residues
	ResidueTypeOPs residue_types_;

	/// for handing out
	ResidueTypeCAPs residue_types_const_;

	/// @brief null list of residues when query fails
	//should make this static or something
	ResidueTypeCAPs empty_residue_list_;


	/// @brief map to ResidueType pointers by AA enum
	std::map< AA, ResidueTypeCAPs > aa_map_;

	/// @brief map to ResidueType pointers by 3-letter string name
	std::map< std::string, ResidueTypeCAPs > name3_map_;

	/// @brief map to ResidueType pointers by unique residue id
	std::map< std::string, ResidueTypeCAP > name_map_;

	/// @brief map to ResidueType pointers by unique residue id, for nonconst access
	std::map< std::string, ResidueTypeOP > nonconst_name_map_;

	/// @brief list of AA types defined
	std::list< AA > aas_defined_;

	//XRW_B_T1
	/*
	/// @brief the coarsify translator in case this a coarse residue set
	coarse::TranslatorSetCOP coarsifier_;
	*/
	//XRW_E_T1

	/// @brief the database directory of the generating files ---> allows to use cached dunbrack libs
	const std::string database_directory_;

private:
	// uncopyable
	ResidueTypeSet( ResidueTypeSet const & );
	ResidueTypeSet const & operator = ( ResidueTypeSet const & );
};

} // chemical
} // core



#endif
