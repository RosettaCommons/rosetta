// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/CacheableAtomID_MapVector.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_core_id_CacheableAtomID_MapVector_hh
#define INCLUDED_core_id_CacheableAtomID_MapVector_hh

// unit headers
#include <core/id/CacheableAtomID_MapVector.fwd.hh>

// type headers
#include <core/types.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// project headers
#include <core/id/AtomID_Map.hh>

// numeric headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace id {

/// @brief AtomID_Map< xyzVector< Real > >
class CacheableAtomID_MapVector : public basic::datacache::CacheableData {
public:
	CacheableAtomID_MapVector();
	virtual ~CacheableAtomID_MapVector();
	virtual basic::datacache::CacheableDataOP clone() const;
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & map() { return map_; }
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > const & map() const { return map_; }
private:
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > map_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace id
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_id_CacheableAtomID_MapVector )
#endif // SERIALIZATION


#endif /* INCLUDED_core_util_datacache_CacheableAtomID_MapVector_HH */
