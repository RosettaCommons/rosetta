// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/WriteableCacheableMap.hh
/// @brief
/// @author Justin Porter


#ifndef INCLUDED_basic_datacache_WriteableCacheableMap_hh
#define INCLUDED_basic_datacache_WriteableCacheableMap_hh

// unit headers
#include <basic/datacache/WriteableCacheableMap.fwd.hh>

// package headers
#include <basic/datacache/WriteableCacheableData.hh>


// C++ headers
#include <map>
#include <string>
#include <set>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <basic/datacache/CacheableData.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {


/// @brief Wrapper for a map< string, string >.
class WriteableCacheableMap : public CacheableData
{
	typedef std::map< std::string, std::set< WriteableCacheableDataOP > > DataMap;

public:
	WriteableCacheableMap();

	WriteableCacheableMap( WriteableCacheableMap const & );

	~WriteableCacheableMap() override;

	CacheableDataOP clone() const override;

	virtual DataMap & map();

	virtual const DataMap & map() const;

	virtual void erase( WriteableCacheableDataOP d );

	virtual DataMap::const_iterator begin() const;

	virtual DataMap::const_iterator end() const;

	virtual std::set< WriteableCacheableDataOP >& operator[] ( std::string const & str );

	virtual DataMap::const_iterator find( std::string const& str ) const;

	virtual bool has( WriteableCacheableDataOP data );

	virtual void insert( WriteableCacheableDataOP data );

	WriteableCacheableMapOP shared_from_this();

private:

	DataMap map_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_WriteableCacheableMap )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_WriteableCacheableMap_HH */
