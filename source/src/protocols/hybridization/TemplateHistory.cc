// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   TemplateHistory.cc
/// @brief
/// @author Frank DiMaio

// Unit headers
#include <protocols/hybridization/TemplateHistory.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

static thread_local basic::Tracer TZ( "protocols.hybridization.TemplateHistory" );

TemplateHistory::TemplateHistory( core::pose::Pose &pose ) {
	// initialize to identity mapping
	core::Size nres = pose.total_residue();
	history_.resize( nres , -1 );
}

int
TemplateHistory::get( core::Size resid ){
	if ( resid <= history_.size() ) {
		return history_[resid];
	} else {
		return -1;
	}
}

void
TemplateHistory::setall( int template_id ) {
	for ( core::Size i=1; i<=history_.size(); ++i ) {
		history_[i] = template_id;
	}
}

void
TemplateHistory::set( core::Size start_res, core::Size stop_res, int template_id) {
	runtime_assert( stop_res<=history_.size() );
	for ( core::Size i=start_res; i<=stop_res; ++i ) {
		history_[i] = template_id;
	}
}

}
//}
}
