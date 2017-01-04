// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableString.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_CacheableString_hh
#define INCLUDED_basic_datacache_CacheableString_hh

// unit headers
#include <basic/datacache/CacheableString.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// C++ headers
#include <string>

#include <platform/types.hh>
#include <utility>
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
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {


/// @brief Wrapper for std::string
class CacheableString : public CacheableData
{
public:
	CacheableString( std::string str );
	~CacheableString() override;
	CacheableDataOP clone() const override;
	virtual std::string const & str() const;

	CacheableStringOP shared_from_this();

private:
	std::string str_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	CacheableString();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_CacheableString )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_CacheableString_HH */
