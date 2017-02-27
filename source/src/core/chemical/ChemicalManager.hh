// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Chemical manager class
///
/// @details
/// The Chemical Manager is a singleton class, which means that it can only been initialized once (exist once in memory). Once initialized,
/// you can call it by simply access it via:
///
/// core::chemical::AtomTypeSetCAP atom_types =
/// core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
///
/// You can substitute AtomTypeSet, with whatever is seen below (residue_type_set, mm_atom_type_set, orbital_type_set).
/// In the below functions, the "tag_in" refers to fullatom, centroid, which basically tells what type of set to load in.
/// The chemical manager will call functions within the AtomTypeSet, MMAtomTypeSet, ResidueTypeSet, etc etc. The classes type set
/// reads in files from the database to create atom types, residue types, and mmatom types. The information from those files are stored
/// in the type class.
///
///
///
/// @author
/// Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// Steven Combs - comments
///
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_ChemicalManager_hh
#define INCLUDED_core_chemical_ChemicalManager_hh

// Unit headers
#include <core/chemical/ChemicalManager.fwd.hh>

// Package headers
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/IdealBondLengthSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/GlobalResidueTypeSet.fwd.hh>

#ifdef MULTI_THREADED

// Utility thread headers
#include <utility/thread/ReadWriteMutex.hh>

#endif

// Utility headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <map>

namespace core {
namespace chemical {

/// @brief a class managing different sets of atom_type_set and residue_type_set
///
/// @details make it as a singleton class so that atom_type_set and residue_type_set are only
///input and initialized once. They can be later retrieved by querying this class.
class ChemicalManager : public utility::SingletonBase< ChemicalManager >
{
public:
	friend class utility::SingletonBase< ChemicalManager >;

public:
	/// @brief query atom_type_set by a name tag
	AtomTypeSetCOP
	atom_type_set( std::string const & tag );

	/// @brief query atom_type_set by a name tag
	ElementSetCOP
	element_set( std::string const & tag );

	//ElementSetCOP
	//gasteiger_element_type_set(std::string const & tag);

	/// @brief query ideal_bond_lengths
	IdealBondLengthSetCOP
	ideal_bond_length_set(std::string const & tag);

	/// @brief query mm_atom_type_set by a name tag
	MMAtomTypeSetCOP
	mm_atom_type_set( std::string const & tag );

	/// @brief query gasteiger_atom_type_set by a name tag
	gasteiger::GasteigerAtomTypeSetCOP
	gasteiger_atom_type_set( std::string const & tag = "default" );

	/// @brief query orbital_type_set by a name tag
	orbitals::OrbitalTypeSetCOP
	orbital_type_set(std::string const & tag);

	/// @brief query residue_type_set by a type
	///
	/// @note If you have access to a Pose/Conformation,
	/// you probably don't want to use this function.
	/// The Conformation can have custom ResidueTypeSets which
	/// add additional ResidueTypes to the ResidueTypeSet.
	/// Instead, use the residue_type_set_for_*() function on the pose.
	/// Those will fall back to this function if there isn't a custom ResidueTypeSet.
	ResidueTypeSetCOP
	residue_type_set( TypeSetMode type_set_type );

	/// @brief query residue_type_set by a name tag
	ResidueTypeSetCOP
	residue_type_set( std::string const & tag );

	/// @brief Check if residue_type_set is instantiated...
	bool has_residue_type_set( std::string const & tag ); //can't be const due to mutex changing

private:
	typedef std::map< std::string, AtomTypeSetOP > AtomTypeSets;
	typedef std::map< std::string, ElementSetOP > ElementSets;
	typedef std::map< std::string, IdealBondLengthSetOP> IdealBondLengthSets;
	typedef std::map< std::string, orbitals::OrbitalTypeSetOP > OrbitalTypeSets;
	typedef std::map< std::string, MMAtomTypeSetOP > MMAtomTypeSets;
	typedef std::map< std::string, gasteiger::GasteigerAtomTypeSetOP > GasteigerAtomTypeSets;
	typedef std::map< std::string, GlobalResidueTypeSetOP > GlobalResidueTypeSets;

#ifdef MULTI_THREADED

private:
	utility::thread::ReadWriteMutex elem_mutex_;
	utility::thread::ReadWriteMutex atomtype_mutex_;
	utility::thread::ReadWriteMutex orbtype_mutex_;
	utility::thread::ReadWriteMutex mmatomtype_mutex_;
	utility::thread::ReadWriteMutex restype_mutex_;
	utility::thread::ReadWriteMutex idealbondlength_mutex_;

#endif

private:

	/// @brief private constructor
	ChemicalManager();

	/// @brief Go and create an atom type set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	AtomTypeSetOP create_atom_type_set( std::string const & tag ) const;

	/// @brief Go and create an element type set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	ElementSetOP create_element_set( std::string const & tag ) const;

	/// @brief Go and create an orbital type set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	orbitals::OrbitalTypeSetOP create_orbital_type_set( std::string const & tag ) const;

	/// @brief Go and create an mm atom type set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	MMAtomTypeSetOP create_mm_atom_type_set( std::string const & tag ) const;

	/// @brief Go and create a residue type set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	GlobalResidueTypeSetOP create_residue_type_set( std::string const & tag ) const;

	/// @brief Go and create an ideal bond length set.  Should be called only after it's been
	/// determined safe (and neccessary) to construct it.
	IdealBondLengthSetOP
	create_ideal_bond_length_set( std::string const & tag ) const;

private: // data
	/// @brief lookup map for querying atom_type_set by name tag
	AtomTypeSets atom_type_sets_;
	/// @brief lookup map for querying element_type_set by name tag
	ElementSets element_sets_;

	/// @brief lookup map for querying orbital_type_set by name tag.
	OrbitalTypeSets orbital_type_sets_;
	/// @brief lookup map for querying mm_atom_type_set by name tag
	MMAtomTypeSets mm_atom_type_sets_;
	/// @brief lookup map for querying gasteiger_atom_type_set by name tag
	GasteigerAtomTypeSets gasteiger_atom_type_sets_;
	/// @brief lookup map for querying residue_type_set by name tag
	GlobalResidueTypeSets residue_type_sets_;

	/// @brief lookup map for the set of ideal bond lengths
	IdealBondLengthSets ideal_bond_length_sets_;
};

/// @details If fail is true, utility_exit if the mode cannot be converted,
/// if not, return the invalid type.
TypeSetMode
type_set_mode_from_string( std::string const & mode, bool fail = true );

std::string
string_from_type_set_mode( TypeSetMode mode );

std::ostream &
operator <<( std::ostream & out, TypeSetMode mode );

} // namespace core
} // namespace chemical

#endif
