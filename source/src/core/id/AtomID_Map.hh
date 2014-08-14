// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/AtomID_Map.hh
/// @brief  Map from Atom identifiers to contained values class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Implemented as a vector< vector > for fast lookup but this has slower
///      insertion/deletion than a std::map or other associative containers
///  @li The outer vector is indexed by the residue number
///  @li The inner vector is indexed by the atom number within the residue
///  @li The map can be sized by first calling resize( n_res ) and then calling
///      resize( i_res, n_atom ) for each residue to set the number of atoms
///  @li When the Value type (T) is bool note that the vector< bool > specialization is
///      used so the values returned by the indexing lookups are not actually references


#ifndef INCLUDED_core_id_AtomID_Map_hh
#define INCLUDED_core_id_AtomID_Map_hh


// Unit headers
#include <core/id/AtomID_Map.fwd.hh>

// Package headers
#include <core/id/AtomID.hh>

#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

#include <utility/vector1_bool.hh>



namespace core {
namespace id {


/// @brief Map from Atom identifiers to contained values class
template< typename T >
class AtomID_Map
{


public: // Types


	typedef  utility::vector1< T >  AtomMap;
	typedef  utility::vector1< AtomMap >  ResidueMap;

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

	/// @brief Default constructor with no arguments (PyRosetta workaround)
	inline
	explicit
	AtomID_Map() :
		default_value_( Value() )
	{}


	/// @brief Default constructor
	inline
	explicit
	AtomID_Map( Value const & default_value_a ) :  // AtomID_Map( Value const & default_value_a = Value() )
		default_value_( default_value_a )
	{}


	/// @brief Number of residues constructor
	inline
	explicit
	AtomID_Map(
		Size const n_res
	) :
		default_value_( Value() ),
		res_map_( n_res )
	{}

	/// @brief Number of residues constructor
	inline
	explicit
	AtomID_Map(
		Size const n_res,
		Value const & default_value_a
	) :
		default_value_( default_value_a ),
		res_map_( n_res )
	{}


	/// @brief Destructor
	inline
	~AtomID_Map()
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
	resize( Size const i_res, Size const n_atom )
	{
		res_map_[ i_res ].resize( n_atom, default_value_ );
	}


	/// @brief Resize the number of atoms of a residue and use a specified fill value
	inline
	void
	resize( Size const i_res, Size const n_atom, Value const & value )
	{
		res_map_[ i_res ].resize( n_atom, value );
	}


	/// @brief Fill the map with the default fill value
	inline
	void
	fill()
	{
		for ( Size i = 1, ie = res_map_.size(); i <= ie; ++i ) {
			AtomMap & atom_map( res_map_[ i ] );
			for ( Size j = 1, je = atom_map.size(); j <= je; ++j ) { // std::fill_n could do this too
				atom_map[ j ] = default_value_;
			}
		}
	}


	/// @brief Fill the map with a specified fill value
	inline
	void
	fill_with( Value const & value )
	{
		for ( Size i = 1, ie = res_map_.size(); i <= ie; ++i ) {
			AtomMap & atom_map( res_map_[ i ] );
			for ( Size j = 1, je = atom_map.size(); j <= je; ++j ) { // std::fill_n could do this too
				atom_map[ j ] = value;
			}
		}
	}


	/// @brief Fill the map at position seqpos with a specified fill value
	inline
	void
	fill_with( Size const seqpos, Value const & value )
	{
		AtomMap & atom_map( res_map_[ seqpos ] );
		for ( Size j = 1, je = atom_map.size(); j <= je; ++j ) { // std::fill_n could do this too
			atom_map[ j ] = value;
		}
	}


	/// @brief Get the value for an AtomID: Extend the map if necessary, filling with the default value
	/// Phil changing this to be a non-resizing function
	inline
	ConstReference
	get( AtomID const & id ) const
	{
		if ( Size( id.rsd() ) > res_map_.size() ) return default_value_;
		//res_map_.resize( id.rsd() );
		AtomMap const & atom_map( res_map_[ id.rsd() ] );
		if ( Size( id.atomno() ) > atom_map.size() ) return default_value_;
		//atom_map.resize( id.atomno(), default_value_ );
		return atom_map[ id.atomno() ];
	}


	/// @brief Set the value for an AtomID: Extend the map if necessary, filling with the default value
	inline
	void
	set( AtomID const & id, Value const & value )
	{
		if ( Size( id.rsd() ) > res_map_.size() ) res_map_.resize( id.rsd() );
		AtomMap & atom_map( res_map_[ id.rsd() ] );
		if ( Size( id.atomno() ) > atom_map.size() ) atom_map.resize( id.atomno(), default_value_ );
		atom_map[ id.atomno() ] = value;
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
		for ( Size i = 1, e = res_map_.size(); i <= e; ++i ) {
			res_map_[ i ].shrink();
		}
		res_map_.shrink();
	}


	/// @brief swap( AtomID_Map )
	inline
	void
	swap( AtomID_Map & s )
	{
		res_map_.swap( s.res_map_ );
	}



	/// @brief swap( AtomID_Map, AtomID_Map )
	template< typename TF >
	friend
	void
	swap( AtomID_Map<TF> & a, AtomID_Map<TF> & b );



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
		AtomID_Map replacement( new_size, default_value_ );
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


	/// @brief Is an element with this AtomID present?
	inline
	bool
	has( AtomID const & id ) const
	{
		return ( ( id.rsd() >= 1 ) && ( id.atomno() >= 1 ) &&
						 ( Size( id.rsd() ) <= res_map_.size() ) && ( Size( id.atomno() ) <= res_map_[ id.rsd() ].size() ) );
	}


public: // Indexers


	// note These do not resize the map


	/// @brief AtomID_Map[ atom_id ] const
	inline
	ConstReference
	operator []( AtomID const & id ) const
	{
		return res_map_[ id.rsd() ][ id.atomno() ];
	}


	/// @brief AtomID_Map[ atom_id ]
	inline
	Reference
	operator []( AtomID const & id )
	{
		return res_map_[ id.rsd() ][ id.atomno() ];
	}


	/// @brief AtomID_Map( atom_id ) const
	inline
	ConstReference
	operator ()( AtomID const & id ) const
	{
		return res_map_[ id.rsd() ][ id.atomno() ];
	}


	/// @brief AtomID_Map( atom_id )
	inline
	Reference
	operator ()( AtomID const & id )
	{
		return res_map_[ id.rsd() ][ id.atomno() ];
	}


	/// @brief AtomID_Map( i_res, i_atom ) const
	inline
	ConstReference
	operator ()( Size const i_res, Size const i_atom ) const
	{
		return res_map_[ i_res ][ i_atom ];
	}


	/// @brief AtomID_Map( i_res, i_atom )
	inline
	Reference
	operator ()( Size const i_res, Size const i_atom )
	{
		return res_map_[ i_res ][ i_atom ];
	}


	/// @brief AtomID_Map[ i_res ] const
	inline
	AtomMap const &
	operator []( Size const i_res ) const
	{
		return res_map_[ i_res ];
	}


	/// @brief AtomID_Map[ i_res ]
	inline
	AtomMap &
	operator []( Size const i_res )
	{
		return res_map_[ i_res ];
	}


	/// @brief AtomID_Map( i_res ) const
	inline
	AtomMap const &
	operator ()( Size const i_res ) const
	{
		return res_map_[ i_res ];
	}


	/// @brief AtomID_Map( i_res )
	inline
	AtomMap &
	operator ()( Size const i_res )
	{
		return res_map_[ i_res ];
	}


public: // Comparison


	/// @brief AtomID_Map == AtomID_Map
	friend
	inline
	bool
	operator ==( AtomID_Map const & a, AtomID_Map const & b )
	{
		return ( a.res_map_ == b.res_map_ );
	}


	/// @brief AtomID_Map != AtomID_Map
	friend
	inline
	bool
	operator !=( AtomID_Map const & a, AtomID_Map const & b )
	{
		return ( a.res_map_ != b.res_map_ );
	}


private: // Fields

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
			ar & default_value_;
			ar & res_map_;
	}
#endif

	/// @brief Default value
	Value default_value_;

	/// @brief Map from Atom identifiers to values
	ResidueMap res_map_;


}; // AtomID_Map


/// @brief swap( AtomID_Map, AtomID_Map )
template< typename T >
inline
void
swap( AtomID_Map<T> & a, AtomID_Map<T> & b )
{
	a.res_map_.swap( b.res_map_ );
}


// PyRosetta concrete types
//class AtomID_Map_bool: public AtomID_Map< bool > {};
class AtomID_Map_Real: public AtomID_Map< Real > {};
class AtomID_Map_AtomID: public AtomID_Map< AtomID > {};


} // namespace id
} // namespace core


#endif // INCLUDED_core_id_AtomID_Map_HH
