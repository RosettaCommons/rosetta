// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin AtomTypeSet
///
/// @brief
/// A class for reading in the atom type properties
///
/// @detailed
/// This class reads in the atom_properties.txt file which contains the "chemical" information for atoms.
/// This does not contain the actual properties, but sets the properties through the AtomType class.
/// This class is called by the ChemicalManager
///
///
///
/// @authors
/// Phil Bradley
/// Steven Combs - comments
///
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_AtomTypeSet_hh
#define INCLUDED_core_chemical_AtomTypeSet_hh


// Unit headers
#include <core/chemical/AtomTypeSet.fwd.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>


// C++ headers
#include <map>

#include <core/types.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

namespace core {
namespace chemical {


/// @brief a set of AtomTypes
///
/// @details a vector of pointers each of which points to an AtomType
/// and the vector index is looked up by an atom_name string in a map
///
class AtomTypeSet : public utility::pointer::ReferenceCount {

public:
	/// @brief c-tor from directory in the rosetta_database This will go
	/// through the directory, usually
	/// "$ROSETTA3_DB/chemical/atom_type_sets/<atom_type_set_name>" and
	/// initialize all the atom types defined in atom_properties.txt and
	/// the extra parameters specified in extras.txt
	AtomTypeSet( std::string const & directory );

	// Construct the atom type set from an sql database
	AtomTypeSet(
		std::string const & name,
		utility::sql_database::sessionOP db_session);

	virtual ~AtomTypeSet();

public:

	/// @brief the name of the AtomTypeSet
	std::string
	name() const;

	/// @brief number of atom types in the set
	Size
	n_atomtypes() const
	{
		return atoms_.size();
	}

	/// @brief  Get the source directory, eg to open an additional file in that directory
	std::string const &
	directory() const
	{
		return directory_;
	}

	/// @brief Check if atom is present
	bool
	has_atom( std::string const & atom_type_name ) const
	{
		return ( atom_type_index_.find( atom_type_name ) != atom_type_index_.end() );
	}

	/// @brief lookup the atom_type by the atom_type_name string
	int
	atom_type_index( std::string const & atom_type_name ) const;

	/// @brief [ ] operator, simulating vector index behavior
	///
	/// @details look up an AtomTypeSet by 1-based indexing
	///
	AtomType const &
	operator[] ( Size const index ) const
	{
		return *( atoms_[ index ] );
	}

	/// @brief [ ] operator, simulating vector index behavior, non-const version
	///
	/// @details look up an AtomTypeSet by 1-based indexing
	///
	AtomType &
	operator[] ( Size const index )
	{
		return *( atoms_[ index ] );
	}

	/// SLOW
	int
	extra_parameter_index( std::string const & name ) const
	{
		std::map< std::string, int >::const_iterator iter( extra_parameter_indices_.find( name ) );
		if ( iter == extra_parameter_indices_.end() ) {
			utility_exit_with_message( "AtomTypeSet: unrecognized atom parameter: "+name );
		}
		return iter->second;
	}

	std::map< std::string, int> const &
	extra_parameter_indices() const {
		return extra_parameter_indices_;
	}

	bool
	has_extra_parameter( std::string const & name ) const
	{
		return ( extra_parameter_indices_.find( name ) != extra_parameter_indices_.end() );
	}

	/// @brief file I/O
	void
	read_file( std::string const & filename );

	/// @brief additional file I/O
	void
	add_parameters_from_file( std::string const & filename );

	// data
private:
	/// @brief  Private helper fxn for performing default parameter substitutions while reading a params file
	Real
	get_default_parameter( std::string const & param_name, std::string const & atm_name ) const;


private: // helper methods for creating an atom type set from a database

	///@brief create an atom type instance on the stack
	AtomType &
	create_atom_type_from_database(
		std::string const & atom_type_set_name,
		std::string const & atom_type_name,
		utility::sql_database::sessionOP db_session);


	void
	read_atom_type_properties_table(
		std::string const & atom_type_set_name,
		chemical::AtomType & atom_type,
		utility::sql_database::sessionOP db_session);

	void
	read_atom_type_extra_parameters_table(
		std::string const & atom_type_set_name,
		chemical::AtomType & atom_type,
		utility::sql_database::sessionOP db_session);

	void
	clone_atom_types_from_commandline();

private:
	/// lookup map: get atom_type_index by atom_type_name
	std::map< std::string, int > atom_type_index_;

	/// @brief a collection of AtomTypes,
	///
	/// @details AtomType has data of atom properties, and it can be
	/// looked up by atom_type_index.
	utility::vector1< AtomType* > atoms_;

	/// @brief lookup map: get atom extra parameter index by atom_type_name
	std::map< std::string, int > extra_parameter_indices_;

	/// @brief  Save the directory name for future use, eg to load associated AtomVDW data
	std::string directory_;

};

} // chemical
} // core

#endif
