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

