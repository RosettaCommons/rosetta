// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/AveragePathLengthFilter.hh
/// @brief Filter on the average covalent path length between residues, including disulfides
/// @author Gabriel Rocklin (grocklin@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_AveragePathLengthFilter_hh
#define INCLUDED_protocols_simple_filters_AveragePathLengthFilter_hh

//unit headers
#include <protocols/simple_filters/AveragePathLengthFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

namespace protocols {
namespace simple_filters {

class AveragePathLengthFilter : public filters::Filter
{
public:
	//default ctor
	AveragePathLengthFilter();

	AveragePathLengthFilter(core::Real path_tightness, core::Real max_path_length)

	;

	AveragePathLengthFilter(
		AveragePathLengthFilter const & src
	);

	virtual ~AveragePathLengthFilter();

	bool
	apply(
		core::pose::Pose const & pose
	) const;

	filters::FilterOP
	clone() const;

	filters::FilterOP
	fresh_instance() const;

	void
	report(
		std::ostream & out,
		core::pose::Pose const & pose
	) const;

	core::Real
	report_sm(
		core::pose::Pose const & pose
	) const;

	core::Real
	compute(
		core::pose::Pose const & pose
	) const;

	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const &
	);

	core::Real
	max_path_length() const;

	void
	max_path_length(
		core::Real value
	);

	core::Real
	path_tightness() const;

	void
	path_tightness(
		core::Real value
	);




private:
	core::Real path_tightness_;
	core::Real max_path_length_;


};

}
}

#endif
