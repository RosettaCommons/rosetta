// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeSet.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_ResidueTypeSet_hh
#define INCLUDED_core_chemical_ResidueTypeSet_hh


// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>

// Package headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueSelector.fwd.hh>

// STL headers
#include <list>

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
	typedef	std::map< std::string, ResidueTypeCOP >::const_iterator const_residue_iterator;

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
		orbitals::OrbitalTypeSetCAP orbital_types
	);

	void
	read_files(
		utility::vector1< std::string > const & filenames,
		AtomTypeSetCAP atom_types,
		ElementSetCAP elements,
		MMAtomTypeSetCAP mm_atom_types,
		orbitals::OrbitalTypeSetCAP orbital_types
	);


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

	void
	add_residue_type(std::string const & tag, std::string const &  filename);

	/// @brief delete a residue type from the set (Use with care)
	void
	remove_residue_type(std::string const & name);

	/// @brief Create correct combinations of adducts for a residue type
	void create_adduct_combinations(
		ResidueType const & rsd,
		std::map< std::string, int > ref_map,
		std::map< std::string, int > count_map,
		utility::vector1< bool > add_mask,
		utility::vector1< Adduct >::const_iterator work_iter
	);

	/// @brief query ResidueTypes by their interchangeability-group name
	ResidueTypeCOPs const &
	interchangeability_group_map( std::string const & name ) const;

	/// @brief query ResidueTypes by their 3-letter name
	ResidueTypeCOPs const &
	name3_map( std::string const & name ) const;

	/// @brief query ResidueType by its unique residue id.
	///
	/// @details since residue id is unique, it only returns
	/// one residue type or exit without match.
	ResidueType const &
	name_map( std::string const & name ) const;

	// IF YOU COMMENT THIS BACK IN, TELL ME WHY BY EMAIL: aleaverfay@gmail.com
	///// @brief query ResidueType by its unique residue id returning a non-const reference to that ResidueType
	///// (Is this really a good idea?  Who is using this functionality?  ResidueTypes are supposed
	///// to be constant once deposited inside a ResidueTypeSet or lots of downstream assumptions
	///// are invalidated!)
	//ResidueType &
	//nonconst_name_map( std::string const & name );

	/// @brief query if a ResidueType of the unique residue id (name) is present.
	bool has_name( std::string const & name ) const;

	/// @brief query if any ResidueTypes in the set have a "name3" tat matches the input name3
	bool
	has_name3( std::string const & name3 ) const;

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
	ResidueTypeCOPs const &
	aa_map( AA const & aa ) const;

	/// @brief select a set of ResidueTypes give certain criteria
	void
	select_residues(
		ResidueSelector const & selector,
		ResidueTypeCOPs & matches
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
	ResidueTypeCOPs const &
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
	add_residue_type_to_maps( ResidueTypeOP rsd );

	void
	remove_residue_type_from_maps(ResidueTypeOP rsd);

	//////////////////
	// data
private:

	/// What does the ChemicalManager call this ResidueTypeSet?
	std::string name_;

	/// @brief the residues
	ResidueTypeOPs residue_types_;

	/// for handing out
	ResidueTypeCOPs residue_types_const_;

	/// @brief null list of residues when query fails
	//should make this static or something
	ResidueTypeCOPs empty_residue_list_;


	/// @brief map to ResidueType pointers by AA enum
	std::map< AA, ResidueTypeCOPs > aa_map_;

	/// @brief map between ResidueType's interchangeability group string and a vector of ResidueTypeCOPs
	std::map< std::string, ResidueTypeCOPs > interchangeability_group_map_;

	/// @brief map to ResidueType pointers by 3-letter string name
	std::map< std::string, ResidueTypeCOPs > name3_map_;

	/// @brief map to ResidueType pointers by unique residue id
	std::map< std::string, ResidueTypeCOP > name_map_;

	/// @brief map to ResidueType pointers by unique residue id, for nonconst access
	std::map< std::string, ResidueTypeOP > nonconst_name_map_;

	/// @brief list of AA types defined
	std::list< AA > aas_defined_;

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
