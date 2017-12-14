// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/DiagnosticData.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_DiagnosticData_hh
#define INCLUDED_basic_datacache_DiagnosticData_hh

// unit headers
#include <basic/datacache/DiagnosticData.fwd.hh>

// type headers

// package headers
#include <basic/datacache/CacheableData.hh>

// C++ headers
#include <map>
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


/// @brief Wrapper for std::map<string, Real>
class DiagnosticData : public CacheableData
{
public:
	DiagnosticData( std::map < std::string, double > const & data_in );
	~DiagnosticData() override;
	CacheableDataOP clone() const override;
	virtual std::map < std::string, double > const & data() const;

	DiagnosticDataOP shared_from_this();

private:
	std::map < std::string, double > data_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DiagnosticData();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( basic_datacache_DiagnosticData )
#endif // SERIALIZATION


#endif /* INCLUDED_basic_datacache_DiagnosticData_HH */
