// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableStringIntegerMap.hh
/// @brief
/// @author Phil Bradley, Jared Adolf-Bryfogle (StringIntMap)


#ifndef INCLUDED_basic_datacache_CacheableStringIntegerMap_hh
#define INCLUDED_basic_datacache_CacheableStringIntegerMap_hh

// unit headers
#include <basic/datacache/CacheableStringIntegerMap.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// C++ headers
#include <map>
#include <string>

#include <platform/types.hh>
#include <utility/down_cast.hh>
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

/// @brief Wrapper for std::map< std::string, integer >
class CacheableStringIntegerMap : public CacheableData
{
public:
	CacheableStringIntegerMap();

	~CacheableStringIntegerMap() override;

	CacheableDataOP
	clone() const override;

	CacheableStringIntegerMapOP
	shared_from_this();

	virtual std::map< std::string, int > &
	map();

	virtual const std::map< std::string, int > &
	map() const;

private:
	std::map< std::string, int > map_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_CacheableStringIntegerMap )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_CacheableStringIntegerMap_HH */
