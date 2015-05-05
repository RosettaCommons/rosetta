// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableStringMap.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_CacheableStringMap_hh
#define INCLUDED_basic_datacache_CacheableStringMap_hh

// unit headers
#include <basic/datacache/CacheableStringMap.fwd.hh>

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


/// @brief Wrapper for a map< string, string >.
class CacheableStringMap : public CacheableData
{
public:
	CacheableStringMap() : CacheableData() {}

	virtual ~CacheableStringMap() {}

	virtual CacheableDataOP clone() const {
		return CacheableDataOP( new CacheableStringMap(*this) );
	}

	virtual std::map< std::string, std::string > & map() {
		return map_;
	}

	virtual const std::map< std::string, std::string > & map() const {
		return map_;
	}

private:

	std::map< std::string, std::string > map_;
};


} // namespace datacache
} // namespace basic

#endif /* INCLUDED_basic_datacache_CacheableStringMap_HH */
