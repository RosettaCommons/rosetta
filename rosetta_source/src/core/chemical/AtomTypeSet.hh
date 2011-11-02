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
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_AtomTypeSet_hh
#define INCLUDED_core_chemical_AtomTypeSet_hh


// Unit headers
#include <core/chemical/AtomTypeSet.fwd.hh>


// Project headers
// AUTO-REMOVED #include <core/chemical/AtomType.hh>

// Utility headers
// Commented by inclean daemon #include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>


// C++ headers
// Commented by inclean daemon #include <string>
#include <map>

#include <core/types.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <utility/vector1_bool.hh>



namespace core {
namespace chemical {


/// @brief a set of AtomTypes
///
/// @details a vector of pointers each of which points to an AtomType
/// and the vector index is looked up by an atom_name string in a map
///
class AtomTypeSet : public utility::pointer::ReferenceCount {

public:
	/// @brief c-tor from directory
	AtomTypeSet( std::string const & directory );

public:

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


	/// @brief lookup the atom_type by the atom_type_name string
	int
	atom_type_index( std::string const & atom_type_name ) const
	{
		std::map< std::string, int >::const_iterator
			iter( atom_type_index_.find( atom_type_name ) );
		if ( iter == atom_type_index_.end() ) {
			utility_exit_with_message("unrecognized atom_type_name "+atom_type_name);
		}
		return iter->second;
	}


	/// @brief [ ] operator, simulating vector index behavior
	///
	/// @details look up an AtomTypeSet by 1-based indexing
	///
	AtomType const &
	operator[] ( Size const index ) const
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




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 	/// is this atom_type an hbond acceptor?
// 	bool
// 	is_acceptor( int const type ) const
// 	{
// 		return atoms_[ type ]->is_acceptor();
// 	}


// 	/// is this atom_type an hbond acceptor?
// 	bool
// 	is_donor( int const type ) const
// 	{
// 		return atoms_[ type ]->is_donor();
// 	}


// 	/// is this atom_type an hbond acceptor?
// 	bool
// 	is_h2o( int const type ) const
// 	{
// 		return atoms_[ type ]->is_h2o();
// 	}


// 	/// is this atom_type a polar_hydrogen?
// 	bool
// 	is_polar_hydrogen( int const type ) const
// 	{
// 		return atoms_[ type ]->is_polar_hydrogen();
// 	}


// 	/// is this atom_type a hydrogen?
// 	bool
// 	is_hydrogen( int const type ) const
// 	{
// 		return atoms_[ type ]->is_hydrogen();
// 	}


// 	///
// 	Real
// 	lk_lambda( int const type ) const
// 	{
// 		return atoms_[ type ]->lk_lambda();
// 	}


// 	///
// 	Real
// 	lk_dgfree( int const type ) const
// 	{
// 		return atoms_[ type ]->lk_dgfree();
// 	}

// 	///
// 	Real
// 	lk_volume( int const type ) const
// 	{
// 		return atoms_[ type ]->lk_volume();
// 	}

// 	///
// 	Real
// 	lj_radius( int const type ) const
// 	{
// 		return atoms_[ type ]->lj_radius();
// 	}

// 	///
// 	Real
// 	lj_wdepth( int const type ) const
// 	{
// 		return atoms_[ type ]->lj_wdepth();
// 	}
