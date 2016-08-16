// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/BasicDataCache.hh
/// @brief  A DataCache storing objects derived from
///         basic::datacache::CacheableData.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_basic_datacache_BasicDataCache_hh
#define INCLUDED_basic_datacache_BasicDataCache_hh

// unit headers
#include <basic/datacache/BasicDataCache.fwd.hh>

// package headers
#include <basic/datacache/DataCache.hh>  // for DataCache
#ifdef WIN32
#include <basic/datacache/CacheableData.hh>
#endif

// C++ headers
#include <cstddef>                       // for size_t

#ifdef    SERIALIZATION
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace basic {
namespace datacache {

/// @brief A DataCache storing objects derived from
///  basic::datacache::CacheableData.
/// @details See DataCache base class for usage details.
class BasicDataCache : public DataCache< CacheableData > {


private: // typedefs


	typedef DataCache< CacheableData > Super;


public: // typedefs


	// typedef Super::size_t size_t;


public: // construct/destruct


	/// @brief default constructor
	BasicDataCache();


	/// @brief size constructor
	/// @param[in] n_types The number of slots for this DataCache.
	BasicDataCache( std::size_t const n_slots );


	/// @brief copy constructor
	BasicDataCache( BasicDataCache const & rval );


	/// @brief default destructor
	virtual
	~BasicDataCache();



public: // assignment


	/// @brief copy assignment
	BasicDataCache & operator =( BasicDataCache const & rval );

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
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_BasicDataCache )
#endif // SERIALIZATION

#endif /* INCLUDED_basic_datacache_BasicDataCache_HH */
