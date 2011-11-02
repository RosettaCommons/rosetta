// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/filters/ScoreCutoffFilter.hh
/// @brief header file for ScoreCutoffFitler class.
/// @detailed
/// @author Florian Richter floric@u.washington.edu


#ifndef INCLUDED_protocols_filters_ScoreCutoffFilter_hh
#define INCLUDED_protocols_filters_ScoreCutoffFilter_hh

// Unit Headers
#include <protocols/filters/ScoreCutoffFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
#include <core/scoring/ScoreType.hh>


#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace filters {


class ScoreCutoffFilter : public Filter {

public:
	typedef Filter parent;

public:
	/// c-tor and
	ScoreCutoffFilter();

	ScoreCutoffFilter( core::Real cutoff_in );

	FilterOP clone() const {
		return new ScoreCutoffFilter( *this ); }

	FilterOP fresh_instance() const {
		return new ScoreCutoffFilter(); }

	virtual void report( std::ostream & ostr, core::pose::Pose const & pose ) const;

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &  );


	/// @brief Returns true if the given pose passes the filter, false otherwise.
	/// In this case, the test is the result of the following comparison:
	/// sc <= cutoff
	/// Where cutoff is the cutoff set for this filter, and sc is the value of the
	/// ScoreType from the Pose Energies object.
	virtual
	bool apply( core::pose::Pose const & pose ) const;

	void set_cutoff( core::Real cutoff_in ){
		cutoff_ = cutoff_in;
	}

	void set_cutoff( core::pose::Pose const & pose ) {
		cutoff_ = get_score( pose ); }

	core::Real cutoff() const {
		return cutoff_;
	}

	void add_score_type( core::scoring::ScoreType scotype );

	void set_score_type( core::scoring::ScoreType scotype );

	void set_positions( utility::vector1< core::Size > const & positions ){
		positions_ = positions;
	}

	void set_unweighted( bool init ){ unweighted_ = init; }

	core::Real get_score( core::pose::Pose const & pose ) const;


	virtual std::string name() const {
		return "ScoreCutoffFilter";
	}

	void
	output_residue_pair_energies( std::ostream & ostr, core::pose::Pose const & pose ) const;

private:

	core::Real cutoff_;
	bool report_residue_pair_energies_;

	utility::vector1< core::scoring::ScoreType > score_types_;

	utility::vector1< core::Size > positions_;

	bool total_score_;
	bool unweighted_; // option to NOT use scorefunction weights when accumulating score

};

} // filters
} // protocols

#endif
