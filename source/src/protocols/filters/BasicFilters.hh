// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/filters/BasicFilters.hh
/// @brief header file for very simple Filter classes
/// @details
/// @author Florian Richter, floric@u.washington.edu (feb 09 ), Sarel Fleishman sarelf@u.washington.edu, Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_filters_BasicFilters_hh
#define INCLUDED_protocols_filters_BasicFilters_hh

// Unit Headers
#include <protocols/filters/BasicFilters.fwd.hh>

// Package headers
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/ResId.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

//// C++ headers
#include <string>

namespace protocols {
namespace filters {

class TrueFilter : public Filter {
public:
	TrueFilter() : Filter( "TrueFilter" ) {}
	bool apply( core::pose::Pose const & ) const { return true; }
	FilterOP clone() const { return FilterOP( new TrueFilter ); }
	FilterOP fresh_instance() const { return FilterOP( new TrueFilter ); }
};

class FalseFilter : public Filter {
public:
	FalseFilter() : Filter( "FalseFilter" ) {}
	bool apply( core::pose::Pose const & ) const { return false; }
	FilterOP clone() const { return FilterOP( new FalseFilter ); }
	FilterOP fresh_instance() const { return FilterOP( new FalseFilter ); }
};

class StochasticFilter : public Filter {

public:
	StochasticFilter();
	virtual ~StochasticFilter();
	StochasticFilter( core::Real const confidence );
	bool apply( core::pose::Pose const & pose ) const;
	FilterOP clone() const;
	FilterOP fresh_instance() const;
	void report( std::ostream &, core::pose::Pose const & ) const {}

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & );

private:
	core::Real confidence_;
};

/// @brief Used to define a compound logical statement involving other filters with
/// AND, OR and XOR
class CompoundFilter : public Filter, public protocols::moves::ResId
{
public:
	typedef std::vector< std::pair< FilterOP, boolean_operations > > CompoundStatement;
	typedef CompoundStatement::iterator iterator;
	typedef CompoundStatement::const_iterator const_iterator;

public:
	CompoundFilter();
	virtual ~CompoundFilter();
	CompoundFilter( CompoundStatement const & );
	bool apply( core::pose::Pose const & ) const;
	FilterOP clone() const;
	FilterOP fresh_instance() const;
	void report( std::ostream &, core::pose::Pose const & ) const;
	core::Real report_sm( core::pose::Pose const & ) const;
	bool compute( core::pose::Pose const & ) const;
	void clear();
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	virtual void set_resid( core::Size const resid );
	void invert( bool const inv );
	void set_reset_filters( utility::vector1<FilterOP> const & reset_filters );
	void reset_filters();
	void clear_reset_filters();

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & );

private:
	core::Real threashold_;
	CompoundStatement compound_statement_;
	bool invert_;
	utility::vector1<FilterOP> reset_filters_;
};

/// @brief Used to combine multiple seperate filters into a single filter value
class CombinedFilter : public Filter
{
public:
	typedef std::pair< FilterOP, core::Real > FilterWeightPair;
	typedef utility::vector1< FilterWeightPair > FilterList;

	CombinedFilter();
	virtual ~CombinedFilter();
	bool apply( core::pose::Pose const & ) const;
	FilterOP clone() const;
	FilterOP fresh_instance() const;
	void report( std::ostream &, core::pose::Pose const & ) const;
	core::Real report_sm( core::pose::Pose const & ) const;
	core::Real compute( core::pose::Pose const & ) const;
	void set_reset_filters( utility::vector1<FilterOP> const & reset_filters );
	void reset_filters();
	void clear_reset_filters();

	/// @brief Set the overall filter threshold.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	inline void set_threshold( core::Real const &val ) { threshold_=val; return; };

	/// @brief Add a filter/weight pair to the list of filters that this filter combines.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	inline void add_filter(
		FilterOP filter,
		core::Real const &weight,
		bool const clone_filter=true
	) {
		filterlist_.push_back( FilterWeightPair( (clone_filter ? filter->clone() : filter), weight) );
		return;
	}

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & );

private:
	core::Real threshold_;
	FilterList filterlist_;
	utility::vector1<FilterOP> reset_filters_;
};


/// @brief Apply a sub-mover prior to calculating a filter value
class MoveBeforeFilter : public Filter
{
public:
	MoveBeforeFilter();
	MoveBeforeFilter(moves::MoverOP mover, FilterCOP filter);
	virtual ~MoveBeforeFilter();
	bool apply( core::pose::Pose const & ) const;
	FilterOP clone() const;
	FilterOP fresh_instance() const;
	void report( std::ostream &, core::pose::Pose const & ) const;
	core::Real report_sm( core::pose::Pose const & ) const;
	//No compute(), as it passes everything on to the sub-mover

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & );

private:
	FilterCOP subfilter_;
	moves::MoverOP submover_;
};

/// @brief Evaluate to a value contingent on the evaluation of another filter.
class IfThenFilter : public Filter
{
public:
	IfThenFilter();
	//IfThenFilter(moves::MoverOP mover, FilterCOP filter);
	virtual ~IfThenFilter();
	bool apply( core::pose::Pose const & ) const;
	FilterOP clone() const;
	FilterOP fresh_instance() const;
	void report( std::ostream &, core::pose::Pose const & ) const;
	core::Real report_sm( core::pose::Pose const & ) const;
	core::Real compute( core::pose::Pose const & ) const;

	void threshold( core::Real threshold ) { threshold_ = threshold; }
	/// @brief Set if threshold is an upper (false/default) or lower (true) limit
	void set_lower_threshold( bool floor = false ) { floor_ = floor; }

	/// @brief Add a condition to the test.
	/// If testfilter evaluates true, then this filter evaluates to valuefilter.
	/// If valuefilter is NULL, then return value instead.
	/// Conditions are evaluated in the order they were added.
	void add_condition( FilterCOP testfilter, FilterCOP valuefilter, core::Real value = 0, bool invert = false, core::Real weight = 1 );

	/// @brief Add evaluation if no conditions trigger
	/// If elsefilter is Null, use absolute value value instead.
	void set_else( FilterCOP elsefilter, core::Real value = 0, core::Real elseweight = 1 );

	void parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & );

private:
	utility::vector1< FilterCOP > iffilters_;
	utility::vector1< FilterCOP > thenfilters_;
	utility::vector1< core::Real > values_;
	utility::vector1< core::Real > weights_;
	/// @brief If true, invert the sense of the iffilter test.
	utility::vector1< bool > invert_;

	FilterCOP elsefilter_;
	core::Real elsevalue_;
	core::Real elseweight_;

	core::Real threshold_;
	/// @brief If true, threshold_ is a lower limit, rather than upper limit
	bool floor_;
};

} // filters
} // protocols

#endif
