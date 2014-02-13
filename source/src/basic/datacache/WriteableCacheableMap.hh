// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

namespace basic {
namespace datacache {


///@brief Wrapper for a map< string, string >.
class WriteableCacheableMap : public CacheableData
{
  typedef std::map< std::string, std::set< WriteableCacheableDataOP > > DataMap;

public:
  WriteableCacheableMap() : CacheableData() {}

  WriteableCacheableMap( WriteableCacheableMap const& other ):
    map_( other.map_ )
  {}

  virtual ~WriteableCacheableMap() {}

  virtual CacheableDataOP clone() const {
    return new WriteableCacheableMap(*this);
  }

  virtual DataMap & map() {
    return map_;
  }

  virtual const DataMap & map() const {
    return map_;
  }

  virtual DataMap::const_iterator begin() const {
    return map_.begin();
  }

  virtual DataMap::const_iterator end() const {
    return map_.end();
  }

  virtual std::set< WriteableCacheableDataOP >& operator[]( std::string const& str ){
    return map_[ str ];
  }

  virtual DataMap::const_iterator find( std::string const& str ) const {
    return map_.find( str );
  }

  virtual void insert( WriteableCacheableDataOP data ){
    map_[ data->datatype() ].insert( data );
  }

private:

  DataMap map_;
};


} // namespace datacache
} // namespace basic

#endif /* INCLUDED_basic_datacache_WriteableCacheableMap_HH */
