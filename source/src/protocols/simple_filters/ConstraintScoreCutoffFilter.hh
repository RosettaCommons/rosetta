// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_filters/ConstraintScoreCutoffFilter.hh
/// @brief header file for ConstraintScoreCutoffFitler class.
/// @details
/// @author Florian Richter floric@u.washington.edu


#ifndef INCLUDED_protocols_simple_filters_ConstraintScoreCutoffFilter_hh
#define INCLUDED_protocols_simple_filters_ConstraintScoreCutoffFilter_hh

// Unit Headers
#include <protocols/simple_filters/ConstraintScoreCutoffFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/constraints/Constraint.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>


#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace simple_filters {


class ConstraintScoreCutoffFilter : public protocols::filters::Filter {

public:
	typedef protocols::filters::Filter parent;

public:
	/// c-tor and
	ConstraintScoreCutoffFilter();
	ConstraintScoreCutoffFilter( core::Real cutoff_in );

	filters::FilterOP clone() const {
		return filters::FilterOP( new ConstraintScoreCutoffFilter( *this ) ); }

	filters::FilterOP fresh_instance() const {
		return filters::FilterOP( new ConstraintScoreCutoffFilter() ); }

	virtual void report( std::ostream & ostr, core::pose::Pose const & pose ) const;
	void parse_my_tag(
										utility::tag::TagCOP tag,
										basic::datacache::DataMap &,
										filters::Filters_map const &,
										protocols::moves::Movers_map const &,
										core::pose::Pose const &
	);

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
		cutoff_ = get_score( pose );
	}

	core::Real cutoff() const {
		return cutoff_;
	}

	void set_constraints( core::scoring::constraints::ConstraintCOPs cst_in );
	void set_score_type( core::scoring::ScoreType scotype );

	core::Real get_score( core::pose::Pose const & pose ) const;
	// Undefined commenting out to fix PyRosetta build  void apply_cst( core::pose::Pose const& pose ) const;

	virtual std::string name() const {
		return "ConstraintScoreCutoffFilter";
	}

private:
	core::scoring::ScoreType score_type_;
	core::Real cutoff_;
	core::scoring::constraints::ConstraintCOPs constraints_;
};

} // filters
} // protocols

#endif
