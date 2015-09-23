// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ScoreTypeFilter.hh
/// @brief definition of filter class ScoreTypeFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Jacob Corn (jecorn@u.washington.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Reworked a bit to make it easier to call this from code, outside of RosettaScripts.

#ifndef INCLUDED_protocols_simple_filters_ScoreTypeFilter_hh
#define INCLUDED_protocols_simple_filters_ScoreTypeFilter_hh

//unit headers
#include <protocols/simple_filters/ScoreTypeFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>

//Auto Headers


namespace protocols {
namespace simple_filters {

class ScoreTypeFilter : public filters::Filter
{
public:
	/// @brief Constructor
	///
	ScoreTypeFilter();

	/// @brief Copy constructor
	///
	ScoreTypeFilter( ScoreTypeFilter const &src );

	/// @brief Constructor with parameters
	///
	ScoreTypeFilter( core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold );

	bool apply( core::pose::Pose const & pose ) const;

	filters::FilterOP clone() const {
		return filters::FilterOP( new ScoreTypeFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const{
		return filters::FilterOP( new ScoreTypeFilter() );
	}

	/// @brief Sets the scorefunction.
	/// @details Note that this filter stores a clone of the scorefunction.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	inline void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in ) { scorefxn_ = sfxn_in->clone(); return; }

	/// @brief Sets the score term that will be used for scoring.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	inline void set_score_type( core::scoring::ScoreType const scoretype_in) { score_type_ = scoretype_in; return; };

	/// @brief Sets the energy threshold used for filtering.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	inline void set_threshold( core::Real const &thresh_in) { score_type_threshold_ = thresh_in; return; }

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	virtual ~ScoreTypeFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	core::Real score_type_threshold_;
	core::scoring::ScoreType score_type_;
	core::scoring::ScoreFunctionOP scorefxn_;
};

}
}

#endif
