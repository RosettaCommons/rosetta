// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/DOF_ID_Map.hh
/// @brief  Map from Atom DOF identifiers to contained values class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Implemented as a vector< AtomID_Map > for fast lookup but this has slower
///      insertion/deletion than a std::map or other associative containers
///  @li The outer vector is indexed by the DOF_Type
///  @li The inner map is indexed by the AtomID (pair of residue number and atom number)
///  @li The DOF layer is on the outside for simple global setting for a DOF and to
///      reuse the AtomID_Map code
///  @li The map can be sized by first calling resize( n_res ) and then calling
///      resize( i_res, n_atom ) for each residue to set the number of atoms
///  @li DOF-specific and uniform default values can be set with the overloaded
///      default_value() functions
///  @li Invariant: The DOF_Map is always of size n_DOF_Type
///  @li When the Value type (T) is bool note that the vector< bool > specialization is
///      used so the values returned by the indexing lookups are not actually references


#ifndef INCLUDED_core_id_DOF_ID_Map_hh
#define INCLUDED_core_id_DOF_ID_Map_hh


// Unit headers
#include <core/id/DOF_ID_Map.fwd.hh>

// Package headers
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>

// Utility headers

#include <core/id/AtomID_Map.fwd.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace id {


/// @brief Map from Atom DOF identifiers to contained values class
template< typename T >
class DOF_ID_Map
{


public: // Types


	typedef  AtomID_Map< T >  AtomMap;
	typedef  utility::vector1< AtomMap >  DOF_Map;

	// STL/boost style
	typedef  T  value_type;
	typedef  typename AtomMap::reference  reference;
	typedef  typename AtomMap::const_reference  const_reference;
	typedef  typename AtomMap::size_type  size_type;

	// Project style
	typedef  T  Value;
	typedef  typename AtomMap::Reference  Reference;
	typedef  typename AtomMap::ConstReference  ConstReference;
	typedef platform::Size Size;


public: // Creation


	/// @brief Default constructor
	inline
	explicit
	DOF_ID_Map( Value const & default_value_a = Value() ) :
		dof_map_( n_DOF_Type, AtomMap( default_value_a ) )
	{}


	/// @brief Number of residues constructor
	inline
	explicit
	DOF_ID_Map(
		Size const n_res,
		Value const & default_value_a = Value()
	) :
		dof_map_( n_DOF_Type, AtomMap( n_res, default_value_a ) )
	{}


	/// @brief Destructor
	inline
	~DOF_ID_Map()
	{}


public: // Methods


	/// @brief Resize to a given number of residues
	inline
	void
	resize( Size const n_res )
	{
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].resize( n_res );
		}
	}


	/// @brief Resize the number of atoms of a residue and use the default fill value
	inline
	void
	resize( Size const i_res, Size const n_atom )
	{
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].resize( i_res, n_atom );
		}
	}


	/// @brief Resize the number of atoms of a residue and use a specified fill value
	inline
	void
	resize( Size const i_res, Size const n_atom, Value const & value )
	{
	debug_assert( dof_map_.size() == n_DOF_Type );
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].resize( i_res, n_atom, value );
		}
	}


	/// @brief Finalize after sizing all the vectors
	inline
	void
	finalize()
	{
		shrink();
	}


	/// @brief Shrink the vectors to remove unused capacity
	inline
	void
	shrink()
	{
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].shrink();
		}
	}


	/// @brief swap( DOF_ID_Map )
	inline
	void
	swap( DOF_ID_Map & s )
	{
		dof_map_.swap( s.dof_map_ );
	}


	/// @brief swap( DOF_ID_Map, DOF_ID_Map )
	friend
	inline
	void
	swap( DOF_ID_Map & a, DOF_ID_Map & b )
	{
		a.dof_map_.swap( b.dof_map_ );
	}


	/// @brief Clear the vectors
	inline
	void
	clear()
	{
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].clear();
		}
	}


	/// @brief Set a DOF for all atoms to the default value for that DOF
	inline
	void
	set( DOF_Type const dof )
	{
		dof_map_[ dof ].fill();
	}


	/// @brief Set a DOF for all atoms to a specified value
	inline
	void
	set( DOF_Type const dof, Value const & value )
	{
		dof_map_[ dof ].fill_with( value );
	}


public: // Properties


	/// @brief Size (Number of DOF types)
	inline
	Size
	size() const
	{
		return dof_map_.size();
	}


	/// @brief Number of residues
	inline
	Size
	n_residue() const
	{
		return dof_map_[ 1 ].n_residue();
	}


	/// @brief Number of atoms in a residue
	inline
	Size
	n_atom( Size const i_res ) const
	{
		return dof_map_[ 1 ].n_atom( i_res );
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return dof_map_.empty();
	}


	/// @brief Default value for a DOF
	inline
	Value const &
	default_value( DOF_Type const dof ) const
	{
		return dof_map_[ dof ].default_value();
	}


	/// @brief Set default value for a DOF
	inline
	void
	default_value( DOF_Type const dof, Value const & default_value_a )
	{
		dof_map_[ dof ].default_value( default_value_a );
	}


	/// @brief Set a uniform default value
	inline
	void
	default_value( Value const & default_value_a )
	{
		for ( Size i = 1; i <= n_DOF_Type; ++i ) {
			dof_map_[ i ].default_value( default_value_a );
		}
	}


public: // Indexers


	/// @brief DOF_ID_Map[ dof_id ] const
	inline
	ConstReference
	operator []( DOF_ID const & id ) const
	{
		return dof_map_[ id.type() ][ id.atom_id() ];
	}


	/// @brief DOF_ID_Map[ dof_id ]
	inline
	Reference
	operator []( DOF_ID const & id )
	{
		return dof_map_[ id.type() ][ id.atom_id() ];
	}


	/// @brief DOF_ID_Map( dof_id ) const
	inline
	ConstReference
	operator ()( DOF_ID const & id ) const
	{
		return dof_map_[ id.type() ][ id.atom_id() ];
	}


	/// @brief DOF_ID_Map( dof_id )
	inline
	Reference
	operator ()( DOF_ID const & id )
	{
		return dof_map_[ id.type() ][ id.atom_id() ];
	}


	/// @brief DOF_ID_Map( i_res, i_atom ) const
	inline
	ConstReference
	operator ()( Size const i_res, Size const i_atom, DOF_Type const dof ) const
	{
		return dof_map_[ dof ][ i_res ][ i_atom ];
	}


	/// @brief DOF_ID_Map( i_res, i_atom )
	inline
	Reference
	operator ()( Size const i_res, Size const i_atom, DOF_Type const dof )
	{
		return dof_map_[ dof ][ i_res ][ i_atom ];
	}


	/// @brief DOF_ID_Map[ dof ] const
	inline
	AtomMap const &
	operator []( DOF_Type const dof ) const
	{
		return dof_map_[ dof ];
	}


	/// @brief DOF_ID_Map[ dof ]
	inline
	AtomMap &
	operator []( DOF_Type const dof )
	{
		return dof_map_[ dof ];
	}


	/// @brief DOF_ID_Map( dof ) const
	inline
	AtomMap const &
	operator ()( DOF_Type const dof ) const
	{
		return dof_map_[ dof ];
	}


	/// @brief DOF_ID_Map( dof )
	inline
	AtomMap &
	operator ()( DOF_Type const dof )
	{
		return dof_map_[ dof ];
	}


public: // Comparison


	/// @brief DOF_ID_Map == DOF_ID_Map
	friend
	inline
	bool
	operator ==( DOF_ID_Map const & a, DOF_ID_Map const & b )
	{
		return ( a.dof_map_ == b.dof_map_ );
	}


	/// @brief DOF_ID_Map != DOF_ID_Map
	friend
	inline
	bool
	operator !=( DOF_ID_Map const & a, DOF_ID_Map const & b )
	{
		return ( a.dof_map_ != b.dof_map_ );
	}


private: // Fields


	/// @brief Map from Atom DOF identifiers to values
	DOF_Map dof_map_;


}; // DOF_ID_Map


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_DOF_ID_Map_HH
