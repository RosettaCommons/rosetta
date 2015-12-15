// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/pack/rotamer_set/WaterPackingInfo.hh
/// @brief
/// @author


#ifndef INCLUDED_core_pack_rotamer_set_WaterPackingInfo_hh
#define INCLUDED_core_pack_rotamer_set_WaterPackingInfo_hh


#include <basic/datacache/CacheableData.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.fwd.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.fwd.hh>

// utility headers

// C++

#include <utility/vector1_bool.hh>

#ifdef WIN32
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#endif


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

class WaterPackingInfo : public basic::datacache::CacheableData {

public:

	WaterPackingInfo();

	WaterPackingInfo( WaterPackingInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const;

	WaterAnchorInfo &
	operator[] ( Size const seqpos );

	WaterAnchorInfo const &
	operator[] ( Size const seqpos ) const;


	void
	clear();

private:
	utility::vector1< WaterAnchorInfoOP > data_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_WaterPackingInfo )
#endif // SERIALIZATION


#endif
