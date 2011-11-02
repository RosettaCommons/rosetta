// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/filters/PDDFScoreFilters.hh
/// @brief header file for PDDFScoreFilter class.
/// @detailed Filters poses based on their fit to the experimental (or fake) P(r) data
/// @author Dominik Gront

#ifndef INCLUDED_protocols_filters_PDDFScoreFilter_hh
#define INCLUDED_protocols_filters_PDDFScoreFilter_hh

#include <protocols/filters/Filter.hh>
#include <protocols/filters/PDDFScoreFilter.fwd.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergy.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace filters {

class PDDFScoreFilter : public Filter {
public:
	/// c-tor and d-tor
	PDDFScoreFilter();

	virtual ~PDDFScoreFilter() {}

	FilterOP clone() const {
		return new PDDFScoreFilter( *this ); }

	FilterOP fresh_instance() const{
		return new PDDFScoreFilter();
	}


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	virtual
	bool apply( core::pose::Pose const & pose ) const;

	virtual std::string name() const {
		return "PDDFScoreFilter";
	}

	core::Real cutoff() const {
		return cutoff_;
	}

	void cutoff(core::Real cutoff_value) { cutoff_ = cutoff_value; }

	core::Real recent_score() { return score_value_; }
private:
	protocols::scoring::methods::saxs::PDDFEnergyOP score_;
	mutable core::Real score_value_;
	core::Real cutoff_;
};

} // filters
} // protocols

#endif
