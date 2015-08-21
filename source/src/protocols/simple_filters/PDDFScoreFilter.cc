// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PDDFScoreFilter.cc
/// @brief runs reject or accept filters on pose
/// @details
///   Contains currently: PDDFScoreFilter
///
///
/// @author Dominik Gront

// Unit Headers
#include <protocols/simple_filters/PDDFScoreFilter.hh>

// Project Headers
#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/filters.OptionKeys.gen.hh>

#include <protocols/scoring/methods/saxs/PDDFEnergy.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//// C++ headers
static thread_local basic::Tracer tr( "protocols.simple_filters.PDDFScoreFilter" );

namespace protocols {
namespace simple_filters {

PDDFScoreFilter::PDDFScoreFilter() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	score_ = protocols::scoring::methods::saxs::PDDFEnergyOP( new protocols::scoring::methods::saxs::PDDFEnergy() );
	cutoff_ = basic::options::option[ basic::options::OptionKeys::filters::set_pddf_filter ]();
	score_value_ = cutoff_+1;
}


bool PDDFScoreFilter::apply( core::pose::Pose const & pose ) const {


	score_value_ = score_->evaluate_pddf_energy( pose );
	if ( score_value_ > cutoff_ ) {
		return false;
	}
	return true;
}

} // simple_filters
} // protocols
