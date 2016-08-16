// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/RelativePoseFilter.hh
/// @brief Computes a filter's value on a pose that is modified from one that is read from disk. Useful for computing values for the same sequence across many structures.
/// @author Sarel Fleishman and Shira Warszawski (shiraw1@weizmann.ac.il)

#ifndef INCLUDED_protocols_simple_filters_MultipleSigmoidsFilter_hh
#define INCLUDED_protocols_simple_filters_MultipleSigmoidsFilter_hh

//unit headers
#include <protocols/simple_filters/MultipleSigmoidsFilter.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <limits>
#include <protocols/simple_filters/RelativePoseFilter.fwd.hh>
#include <protocols/simple_filters/OperatorFilter.fwd.hh>
#include <protocols/simple_filters/SigmoidFilter.fwd.hh>


namespace protocols {
namespace simple_filters {

/// @brief simply takes a list of pdbs and creates relative pose then extract sigmoids and call operator (product)
class MultipleSigmoids : public filters::Filter
{
public:
	MultipleSigmoids();
	virtual ~MultipleSigmoids();
	filters::FilterOP clone() const;
	filters::FilterOP fresh_instance() const;
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & );
	core::Real compute( core::pose::Pose const & pose ) const;
	core::Real threshold() const{ return threshold_; }
	void threshold( core::Real const t ){ threshold_ = t; }
	OperatorOP operator_filter() const;
	void operator_filter( OperatorOP opt );
	SigmoidOP sigmoid_filter() const;
	// Undefined, commenting out to fix PyRosetta build  void sigmoid_filter( SigmoidOP sig );
	RelativePoseFilterOP relative_pose_filter() const;
	void relative_pose_filter( RelativePoseFilterOP rpose);
	void reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ); /// allows within-trajectory resetting of the baseline. Notice this is nonconst, so can't be called from apply. attempt_read_from_checkpoint should be true for MC trials > 1, but false otherwise

private:
	std::string file_names_; // dflt ""
	core::Real threshold_; // dflt 0
	RelativePoseFilterOP r_pose_; //dflt NULL
	SigmoidOP sig_; //dflt NULL
	OperatorOP operatorF_; //dflt NULL
};
}
}

#endif
