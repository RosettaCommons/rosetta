// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/NamedAtomID_Map.hh
/// @brief  Map from named Atom identifiers to contained values class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com), Rocco Moretti (rmorettiase@gmail.com)
///
/// @note
///  @li Implemented as a vector< std::map > for fast residue lookup but this has slower
///      insertion/deletion than a std::map or other associative containers
///  @li The outer vector is indexed by the residue number
///  @li The inner map is a map indexed by name strings (white-space stripped)
///  @li The interface tries to conform to that of AtomID_Map, aside from the nessecary atom index -> string conversion
///  @li The two are not related by inheritance, however, just as AtomID and NamedAtomID aren't related.


#ifndef INCLUDED_core_id_NamedAtomID_Map_hh
#define INCLUDED_core_id_NamedAtomID_Map_hh

// Unit headers
#include <core/id/NamedAtomID_Map.fwd.hh>

// Package headers
#include <core/id/NamedAtomID.hh>

// Utility headers

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <map>
#include <string>

namespace core {
namespace id {


/// @brief Map from Atom identifiers to contained values class
template< typename T >
class NamedAtomID_Map
{


public: // Types

	typedef  std::map< std::string, T >  AtomMap;
	typedef  utility::vector1< AtomMap >  ResidueMap;

	// STL/boost style
	typedef  T  value_type;
	typedef  typename AtomMap::mapped_type  reference;
	typedef  typename AtomMap::mapped_type const  const_reference;
	typedef  typename AtomMap::size_type  size_type;

	// Project style
	typedef  T  Value;
	typedef  typename AtomMap::mapped_type  Reference;
	typedef  typename AtomMap::mapped_type const ConstReference;
	typedef platform::Size Size;


public: // Creation

	/// @brief Default constructor with no arguments (PyRosetta workaround)
	inline
	explicit
	NamedAtomID_Map() :
		default_value_( Value() )
	{}


	/// @brief Default constructor
	inline
	explicit
	NamedAtomID_Map( Value const & default_value_a ) :
		default_value_( default_value_a )
	{}


	/// @brief Number of residues constructor
	inline
	explicit
	NamedAtomID_Map(
		Size const n_res
	) :
		default_value_( Value() ),
		res_map_( n_res )
	{}

	/// @brief Number of residues constructor
	inline
	explicit
	NamedAtomID_Map(
		Size const n_res,
		Value const & default_value_a
	) :
		default_value_( default_value_a ),
		res_map_( n_res )
	{}


	/// @brief Destructor
	inline
	~NamedAtomID_Map()
	{}


public: // Methods


	/// @brief Resize to a given number of residues
	inline
	void
	resize( Size const n_res )
	{
		res_map_.resize( n_res );
	}


	/// @brief Resize the number of atoms of a residue and use the default fill value
	inline
	void
	resize( Size const, Size const )
	{
		// No-op - maps automatically resize. Kept for interface compatibility.
	}


	/// @brief Resize the number of atoms of a residue and use a specified fill value
	inline
	void
	resize( Size const, Size const, Value const & )
	{
		// No-op - maps automatically resize. Kept for interface compatibility.
	}


	/// @brief Fill the map with the default fill value
	inline
	void
	fill()
	{
		// Actually just clears entries empty map means a default value.
		for ( Size i = 1, ie = res_map_.size(); i <= ie; ++i ) {
			res_map_[ i ].clear();
		}
	}


	// brief Fill the map with a specified fill value
	// Not supported - can't fill with arbitrary value
	//inline
	//void
	//fill_with( Value const & value ) {}

	// brief Fill the map at position seqpos with a specified fill value
	// Not supported - can't fill with arbitrary value
	//inline
	//void
	//fill_with( Size const seqpos, Value const & value ) {}

	/// @brief Get the value for an NamedAtomID: Return default value if not present
	/// Phil changing this to be a non-resizing function
	inline
	ConstReference
	get( NamedAtomID const & id ) const
	{
		if ( Size( id.rsd() ) > res_map_.size() ) return default_value_;
		AtomMap const & atom_map( res_map_[ id.rsd() ] );
		if (  ! atom_map.count( id.atom() ) )  return default_value_;
		return atom_map.find( id.atom() )->second;
	}


	/// @brief Set the value for an NamedAtomID: Extend the map if necessary, filling with the default value
	inline
	void
	set( NamedAtomID const & id, Value const & value )
	{
		if ( Size( id.rsd() ) > res_map_.size() ) res_map_.resize( id.rsd() );
		AtomMap & atom_map( res_map_[ id.rsd() ] );
		atom_map[ id.atom() ] = value;
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
		res_map_.shrink();
	}


	/// @brief swap( NamedAtomID_Map )
	inline
	void
	swap( NamedAtomID_Map & s )
	{
		res_map_.swap( s.res_map_ );
	}


	/// @brief swap( NamedAtomID_Map, NamedAtomID_Map )
	template< typename TF >
	friend
	void
	swap( NamedAtomID_Map<TF> & a, NamedAtomID_Map<TF> & b );


	/// @brief Clear the map
	inline
	void
	clear()
	{
		// default value isn't changed
		res_map_.clear();
	}


	/// @brief Clear the map and set a new default value
	inline
	void
	clear( Value const & default_value_a )
	{
		default_value_ = default_value_a;
		res_map_.clear();
	}


	/// Should move to .cc?
	/// if old2new[pos] == 0 , that position's mapping is lost
	/// if old2new[1...old_size] doesnt cover all of [1...new_size], the missed positions will have res_map_[pos].empty()
	inline
	void
	update_sequence_numbering( Size const new_size, utility::vector1< Size > const & old2new )
	{
		// swap is very slick
		NamedAtomID_Map replacement( new_size, default_value_ );
		for ( Size i=1, i_end = size(); i<= i_end; ++i ) {
			Size const new_pos( old2new[ i ] );
			if ( new_pos ) {
				replacement[ new_pos ].swap( res_map_[ i ] );
			}
		}
		swap( replacement );
	}

public: // Properties


	/// @brief Size
	inline
	Size
	size() const
	{
		return res_map_.size();
	}


	/// @brief Number of residues (size)
	inline
	Size
	n_residue() const
	{
		return res_map_.size();
	}


	/// @brief Number of atoms in a residue
	inline
	Size
	n_atom( Size const i_res ) const
	{
		return res_map_[ i_res ].size();
	}


	/// @brief Empty?
	inline
	bool
	empty() const
	{
		return res_map_.empty();
	}


	/// @brief Default value
	inline
	Value const &
	default_value() const
	{
		return default_value_;
	}


	/// @brief Set default value
	inline
	void
	default_value( Value const & default_value_a )
	{
		default_value_ = default_value_a;
	}


	/// @brief Is an element with this NamedAtomID present?
	inline
	bool
	has( NamedAtomID const & id ) const
	{
		return ( ( id.rsd() >= 1 ) && ( Size( id.rsd() ) <= res_map_.size() ) &&
			( res_map_[ id.rsd() ].count( id.atom() ) ) );
	}

public: // Indexers

	// note These do not resize the map

	/// @brief NamedAtomID_Map[ atom_id ] const
	inline
	ConstReference
	operator []( NamedAtomID const & id ) const
	{
		return (*this)( id );
	}


	/// @brief NamedAtomID_Map[ atom_id ]
	inline
	Reference
	operator []( NamedAtomID const & id )
	{
		return (*this)( id );
	}


	/// @brief NamedAtomID_Map( atom_id ) const
	inline
	ConstReference
	operator ()( NamedAtomID const & id ) const
	{
		return (*this)( id.rsd(), id.atom() );
	}


	/// @brief NamedAtomID_Map( atom_id )
	inline
	Reference
	operator ()( NamedAtomID const & id )
	{
		return (*this)( id.rsd(), id.atom() );
	}


	/// @brief NamedAtomID_Map( i_res, i_atom ) const
	inline
	ConstReference
	operator ()( Size const i_res, std::string const atom ) const
	{
		if ( res_map_[ i_res ].count( atom ) ) {
			return res_map_[ i_res ].find( atom )->second;
		} else {
			return default_value_;
		}
	}


	/// @brief NamedAtomID_Map( i_res, i_atom )
	inline
	Reference
	operator ()( Size const i_res, std::string const atom )
	{
		//Hack until I find out the root cause of this from file_data.cc
		if ( i_res > res_map_.size() ) {
			throw utility::excn::EXCN_RangeError("Residue outside res_map range");
		}
		if ( res_map_[ i_res ].count( atom ) == 0 ) {
			res_map_[ i_res ][ atom ] = default_value_;
		}
		return res_map_[ i_res ][ atom ];
	}


	/// @brief NamedAtomID_Map[ i_res ] const
	inline
	AtomMap const &
	operator []( Size const i_res ) const
	{
		return res_map_[ i_res ];
	}


	/// @brief NamedAtomID_Map[ i_res ]
	inline
	AtomMap &
	operator []( Size const i_res )
	{
		return res_map_[ i_res ];
	}


	/// @brief NamedAtomID_Map( i_res ) const
	inline
	AtomMap const &
	operator ()( Size const i_res ) const
	{
		return res_map_[ i_res ];
	}


	/// @brief NamedAtomID_Map( i_res )
	inline
	AtomMap &
	operator ()( Size const i_res )
	{
		return res_map_[ i_res ];
	}


public: // Comparison


	/// @brief NamedAtomID_Map == NamedAtomID_Map
	friend
	inline
	bool
	operator ==( NamedAtomID_Map const & a, NamedAtomID_Map const & b )
	{
		return ( a.res_map_ == b.res_map_ );
	}


	/// @brief NamedAtomID_Map != NamedAtomID_Map
	friend
	inline
	bool
	operator !=( NamedAtomID_Map const & a, NamedAtomID_Map const & b )
	{
		return ( a.res_map_ != b.res_map_ );
	}

#ifdef    SERIALIZATION
	/// @brief Serialization routine.  In order for this to successfully compile,
	/// the code that is trying to serialize an instance of this class will need
	/// to #include <utility/serialization/serialization.hh>, which should already
	/// happen, and also #include <utility/vector1.srlz.hh>.
	template < class Archive >
	void save( Archive & arc ) const {
		arc( CEREAL_NVP( default_value_ ), CEREAL_NVP( res_map_ ) );
	}

	/// @brief Deserialization routine
	template < class Archive >
	void load( Archive & arc ) {
		arc( default_value_, res_map_ );
	}
#endif // SERIALIZATION

private: // Fields

	/// @brief Default value
	Value default_value_;

	/// @brief Map from Atom identifiers to values
	ResidueMap res_map_;

}; // NamedAtomID_Map


/// @brief swap( NamedAtomID_Map, NamedAtomID_Map )
template< typename T >
inline
void
swap( NamedAtomID_Map<T> & a, NamedAtomID_Map<T> & b )
{
	a.res_map_.swap( b.res_map_ );
}


// PyRosetta concrete types
class NamedAtomID_Map_NamedAtomID: public NamedAtomID_Map< NamedAtomID > {};

} // namespace id
} // namespace core

#endif // INCLUDED_core_id_NamedAtomID_Map_HH
