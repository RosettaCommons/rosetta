// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   CanonicalFragmentHistory.cc
/// @brief
/// @author TJ Brunette

// Unit headers
#include <core/scoring/methods/vall_lookback/VallLookbackData.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace methods {

static basic::Tracer TZ("protocols.simple_moves.VallLookbackData");

VallLookbackData::VallLookbackData( core::pose::Pose &pose ) {
	// initialize to identity mapping
	core::Size nres = pose.total_residue();
	rmsd_history_.resize( nres , -1 );
	res_changed_history_.resize(nres,true);
}

void
VallLookbackData::set_rmsd( core::Size resid, core::Real rmsd) {
	runtime_assert( resid<=rmsd_history_.size() );
	rmsd_history_[resid] = rmsd;
}

core::Real
VallLookbackData::get_rmsd( core::Size resid ) const{
	if ( resid <= rmsd_history_.size() ) {
		return rmsd_history_[resid];
	} else {
		return -1;
	}
}

void
VallLookbackData::set_res_changed( core::Size resid, bool changed) {
	runtime_assert( resid<=res_changed_history_.size() );
	res_changed_history_[resid] = changed;
}

bool
VallLookbackData::get_res_changed( core::Size resid ) const{
	if ( resid <= res_changed_history_.size() ) {
		return res_changed_history_[resid];
	} else {
		return false;
	}
}

}//end methods
}//end scoring
}//end core


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::methods::VallLookbackData::VallLookbackData() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::methods::VallLookbackData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( rmsd_history_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( res_changed_history_ ) ); // utility::vector1<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::methods::VallLookbackData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( rmsd_history_ ); // utility::vector1<core::Real>
	arc( res_changed_history_ ); // utility::vector1<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::methods::VallLookbackData );
CEREAL_REGISTER_TYPE( core::scoring::methods::VallLookbackData )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_methods_vall_lookback_VallLookbackData )
#endif // SERIALIZATION
