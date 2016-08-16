// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DockFilters
/// @brief Filters for docking
/// @author Jeff Gray


#ifndef INCLUDED_protocols_docking_DockFilters_hh
#define INCLUDED_protocols_docking_DockFilters_hh

// Unit Headers
#include <protocols/docking/types.hh>
#include <protocols/docking/DockFilters.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/ScoreCutoffFilter.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers

namespace protocols {
namespace docking {

/// @brief Low-resolution (centroid-mode) filter for docking.
/// Checks (1) at least some contact is being made between docking partners,
///        (2) clashes are limited so partners are not overlapping
///    and (3) constraints, if present, are met
class DockingLowResFilter : public protocols::filters::Filter
{
public:
	DockingLowResFilter();
	DockingLowResFilter( const DockingLowResFilter & init );
	~DockingLowResFilter();
	void set_use_constraints( bool flag, core::Real cutoff=1.0 ); /// @brief add docking constraints
	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose) const;
	protocols::filters::FilterOP clone() const { return protocols::filters::FilterOP( new DockingLowResFilter( *this ) ); }
	protocols::filters::FilterOP fresh_instance() const { return protocols::filters::FilterOP( new DockingLowResFilter() ); }

private:
	bool use_constraints_;            /// @brief boolean to indicate if constraints are used
	core::Real constraint_cutoff_;    /// @brief cutoff value for the constraint score
	protocols::filters::FilterCollectionOP filters_;
};


/// @brief High-resolution (all-atom) filter for docking.
/// Checks (1) total_score beats the cutoff given
///        (2) interface_score must be negative
/// @details /// @rosetta++ had: 1-score filter 2-fa_rep filter 3-interfaceE filter 4-chainbreak filter
/// TTD: add these other filters
class DockingHighResFilter : public protocols::filters::Filter
{
public:
	DockingHighResFilter();
	DockingHighResFilter( const DockingHighResFilter & init );
	~DockingHighResFilter();
	void set_score_margin( core::Real new_score_margin );
	void set_score_cutoff( core::Real new_cutoff ) { scorefilter_->set_cutoff( new_cutoff ); }
	void set_moveable_jumps( DockJumps const & movable_jumps ) { movable_jumps_ = movable_jumps; }
	void set_scorefunction( core::scoring::ScoreFunctionOP const scorefunction );
	bool apply( core::pose::Pose const & pose ) const;
	//core::Real report_interface_score() const { return interface_score_; } // only valid after apply()
	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const { return protocols::filters::FilterOP( new DockingHighResFilter() ); }

private:
	DockJumps movable_jumps_;
	core::Real score_margin_; /// @brief extra margin for passing filters for early in protocol
	core::scoring::ScoreFunctionOP scorefunction_; /// @brief ScoreFunction for evaluating interface_score //defaults to docking scorefxn
	protocols::simple_filters::ScoreCutoffFilterOP scorefilter_; /// @brief filter for total_score
};


} // docking
} // protocols

#endif
