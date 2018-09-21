// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableUint64MathMatrixFloatMap.hh
/// @brief
/// @author Brian Coventry


#ifndef INCLUDED_basic_datacache_CacheableUint64MathMatrixFloatMap_hh
#define INCLUDED_basic_datacache_CacheableUint64MathMatrixFloatMap_hh

// unit headers
#include <basic/datacache/CacheableUint64MathMatrixFloatMap.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>
#include <numeric/MathMatrix.hh>

// C++ headers
#include <unordered_map>
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

/// @brief Wrapper for std::map< uint64_t, MathMatrix<float> >
class CacheableUint64MathMatrixFloatMap : public CacheableData
{
public:
	CacheableUint64MathMatrixFloatMap();

	~CacheableUint64MathMatrixFloatMap() override;

	CacheableDataOP
	clone() const override;

	CacheableUint64MathMatrixFloatMapOP
	shared_from_this();

	virtual std::unordered_map< uint64_t, numeric::MathMatrix<float> > &
	map();

	virtual const std::unordered_map< uint64_t, numeric::MathMatrix<float> > &
	map() const;

private:
	std::unordered_map< uint64_t, numeric::MathMatrix<float> > map_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_CacheableUint64MathMatrixFloatMap )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_CacheableStringFloatMap_HH */
