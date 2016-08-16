// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/filters/SAXSScoreFilters.hh
/// @brief header file for SAXSScoreFilter class.
/// @details Filters poses based on full atom SAXS energy
/// @author Dominik Gront

#ifndef INCLUDED_protocols_simple_filters_SAXSScoreFilter_hh
#define INCLUDED_protocols_simple_filters_SAXSScoreFilter_hh

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/SAXSScoreFilter.fwd.hh>
#include <core/scoring/saxs/SAXSEnergyFA.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

class SAXSScoreFilter : public protocols::filters::Filter {
public:
	/// c-tor and d-tor
	SAXSScoreFilter();

	virtual ~SAXSScoreFilter() {}

	filters::FilterOP clone() const {
		return filters::FilterOP( new SAXSScoreFilter( *this ) ); }

	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new SAXSScoreFilter() );
	}


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	virtual
	bool apply( core::pose::Pose const & pose ) const;

	virtual std::string name() const {
		return "SAXSScoreFilter";
	}

	core::Real cutoff() const {
		return cutoff_;
	}

	void cutoff(core::Real cutoff_value) { cutoff_ = cutoff_value; }

	core::Real recent_score() { return score_value_; }
private:
	core::scoring::saxs::SAXSEnergyFA score;
	mutable core::Real score_value_;
	core::Real cutoff_;
};

} // filters
} // protocols

#endif
