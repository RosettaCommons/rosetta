// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableStringFloatMap.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_CacheableStringFloatMap_hh
#define INCLUDED_basic_datacache_CacheableStringFloatMap_hh

// unit headers
#include <basic/datacache/CacheableStringFloatMap.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// C++ headers
#include <map>
#include <string>

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


/// @brief Wrapper for std::map< std::string, float >
class CacheableStringFloatMap : public CacheableData
{
public:
	CacheableStringFloatMap() : CacheableData() {}
	~CacheableStringFloatMap() override= default;
	CacheableDataOP clone() const override { return CacheableDataOP( new CacheableStringFloatMap(*this) ); }

	CacheableStringFloatMapOP shared_from_this() { return utility::pointer::static_pointer_cast<CacheableStringFloatMap>( CacheableData::shared_from_this() ); }

	virtual std::map< std::string, float > & map(){ return map_; }
	virtual const std::map< std::string, float > & map() const { return map_; }
private:
	std::map< std::string, float > map_;
};


} // namespace datacache
} // namespace basic

#endif /* INCLUDED_basic_datacache_CacheableStringFloatMap_HH */
