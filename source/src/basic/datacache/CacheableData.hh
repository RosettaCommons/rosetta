// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableData.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_CacheableData_hh
#define INCLUDED_basic_datacache_CacheableData_hh

// unit headers
#include <basic/datacache/CacheableData.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

/// @brief base class for data storable within a DataCache
class CacheableData : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< CacheableData >
{
public:
	/// self pointers
	inline CacheableDataCOP get_self_ptr() const { return shared_from_this(); }
	inline CacheableDataOP get_self_ptr() { return shared_from_this(); }
	inline CacheableDataCAP get_self_weak_ptr() const { return CacheableDataCAP( shared_from_this() ); }
	inline CacheableDataAP get_self_weak_ptr() { return CacheableDataAP( shared_from_this() ); }

	virtual ~CacheableData();

	virtual CacheableDataOP clone() const = 0;

#ifdef    SERIALIZATION
	template < class Archive >
	void save( Archive & archive ) const;

	template < class Archive >
	void load( Archive & archive );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_CacheableData )
#endif // SERIALIZATION

#endif /* INCLUDED_basic_datacache_CacheableData_HH */
