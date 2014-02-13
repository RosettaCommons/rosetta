// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/DataCache.hh
/// @brief  General miscellaneous data storage for a Pose.
/// @author Phil Bradley
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


#ifndef INCLUDED_basic_datacache_DataCache_hh
#define INCLUDED_basic_datacache_DataCache_hh

// unit headers
#include <basic/datacache/DataCache.fwd.hh>

#include <basic/datacache/CacheableData.hh>


// type headers

// utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <vector>
#include <basic/datacache/CacheableData.fwd.hh>



namespace basic {
namespace datacache {


/// @brief Indexed storage for objects derived from a ReferenceCountable
///  data type.
/// @details Intended for use as a generic data cache by storing objects
///  derived from a ReferenceCountable data type in a unique slot designated
///  by an integer id (enum, size index, etc.). The DataCache will only store
///  one object per slot/id.  For example, see the PoseDataCache used in
///  core::pose::Pose, which is indexed by the enum basic::pose::datacache:CacheableDataType.
///  Currently when data is set(), it is not cloned -- classes deriving from
///  DataCache should remember to overload set() if they need cloning behavior.
/// @tparam Data Class derived from utility::pointer::ReferenceCount that
///  defines a virtual clone() method.
template< typename Data >
class DataCache : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // typedefs


	typedef utility::pointer::owning_ptr< Data > DataOP;
	typedef utility::pointer::owning_ptr< Data const > DataCOP;
	typedef utility::pointer::access_ptr< Data > DataAP;
	typedef utility::pointer::access_ptr< Data const > DataCAP;


protected: // typedefs


	typedef utility::vector1< DataOP > DataOPs;


public: // constructors


	/// @brief default constructor
	DataCache() :
		Super()
	{}


	/// @brief size constructor
	/// @param[in] n_types The number of slots for this DataCache.
	DataCache( size_t const n_slots ) :
		Super(),
		data_( n_slots, 0 )
	{}


	/// @brief copy constructor
	DataCache( DataCache const & src ) :
		Super( src )
	{
		*this = src;
	}


	/// @brief destructor
	virtual
	~DataCache() {}


public: // assignment


	/// @brief copy assignment
	DataCache &
	operator =( DataCache const & src )
	{
		if ( this != &src ) {
			Super::operator =( src );

			data_.resize( src.data_.size() );
			for ( size_t i = 1; i <= src.data_.size(); ++i ) {
				if ( src.data_[i] ) {
					data_[i] = src.data_[i]->clone();
				} else {
					data_[i] = 0;
				}
			}
		}

		return *this;
	}


public: // state


	/// @brief the number of slots in this cache
	inline
	size_t
	size() const {
		return data_.size();
	}


	/// @brief resize the cache for the given number of slots
	/// @param[in] n_slots The new number of slots.
	inline
	void
	resize( size_t const n_slots )
	{
		data_.resize( n_slots, 0 );
	}


	/// @brief clear all stored data
	inline
	void
	clear()
	{
		data_.assign( data_.size(), 0 );
	}


	/// @brief clear the object in a selected slot
	inline
	void
	clear( size_t const slot )
	{
		data_[ slot ] = 0;
	}


	/// @brief is there an object in the slot?
	inline
	bool
	has( size_t const slot ) const
	{
		return ( data_[ slot ] != 0 );
	}


public: // accessors


	/// @brief get base class reference to the object stored in the slot
	inline
	Data const &
	get( size_t const slot ) const
	{
		assert( data_[ slot ] );
		return *( data_[ slot ] );
	}


	/// @brief get derived class reference to the object stored in the slot
	/// @tparam D  class derived from Data
	template< typename D >
	inline
	D const &
	get( size_t const slot ) const
	{
		return static_cast< D const & >( get( slot ) );
	}


	/// @brief get base class reference to object stored in the slot
	inline
	Data &
	get( size_t const slot )
	{
		assert( data_[ slot ] );
		return *( data_[ slot ] );
	}


	/// @brief get derived class reference to object stored in the slot
	/// @tparam D  class derived from Data
	template< typename D >
	inline
	D &
	get( size_t const slot )
	{
		return static_cast< D & >( get( slot ) );
	}


	/// @brief get base class owning ptr to object stored in the slot
	inline
	DataCOP
	get_const_ptr( size_t const slot ) const
	{
		assert( data_[ slot ] );
		return data_[ slot ];
	}


	/// @brief get derived class owning ptr to object stored in the slot
	/// @tparam D class derived from Data
	template< typename D >
	inline
	utility::pointer::owning_ptr< D const >
	get_const_ptr( size_t const slot ) const
	{
		return utility::pointer::static_pointer_cast< D const >( get_const_ptr( slot ) );
	}


	/// @brief get base class owning ptr to object stored in the slot
	inline
	DataOP
	get_ptr( size_t const slot )
	{
		assert( data_[ slot ] );
		return data_[ slot ];
	}


	/// @brief get derived class owning ptr to object stored in the slot
	/// @tparam D class derived from Data
	template< typename D >
	inline
	utility::pointer::owning_ptr< D >
	get_ptr( size_t const slot )
	{
		return utility::pointer::static_pointer_cast< D >( get_ptr( slot ) );
	}


	/// @brief get base class raw ptr to object stored in the slot
	inline
	Data const *
	get_raw_const_ptr( size_t const slot ) const
	{
		assert( data_[ slot ] );
		return data_[ slot ].get();
	}


	/// @brief get derived class raw ptr to object stored in the slot
	/// @tparam D class derived from Data
	template< typename D >
	inline
	D const *
	get_raw_const_ptr( size_t const slot ) const
	{
		return static_cast< D const * >( get_raw_const_ptr( slot ) );
	}


	/// @brief get base class raw ptr to object stored in the slot
	inline
	Data *
	get_raw_ptr( size_t const slot )
	{
		assert( data_[ slot ] );
		return data_[ slot ].get();
	}


	/// @brief get derived class raw ptr to object stored in the slot
	/// @tparam D class derived from Data
	template< typename D >
	inline
	D *
	get_raw_ptr( size_t const slot )
	{
		return static_cast< D * >( get_raw_ptr( slot ) );
	}


public: // mutators


	/// @brief store data in the given slot
	/// @param[in] The slot to use.
	/// @param[in] observer The data to store -- data is *not* cloned.
	inline
	void
	set( size_t const slot, DataOP new_data )
	{
		// yab: Switching the assignment below to a clone() clearly causes
		// breakage in the code and tests fail.  Someone needs to sit down
		// and figure out if this is the desired behavior.
		data_[ slot ] = new_data;
	}


protected: // access


	/// @brief get the storage vector
	inline
	DataOPs const &
	data() const {
		return data_;
	}


	/// @brief get the storage vector
	inline
	DataOPs &
	data() {
		 return data_;
	}


private: // data


	/// @brief stores data via unique integer index
	DataOPs data_;


};

// PyRosetta workaround
class DataCache_CacheableData : public DataCache<CacheableData> {};


} // namespace datacache
} // namespace basic


#endif /* INCLUDED_basic_datacache_DataCache_HH */
