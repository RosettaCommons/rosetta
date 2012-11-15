// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueCountFilter.hh
/// @brief Filter on the total number of residues in the structure
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_ResidueCountFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueCountFilter_hh

//unit headers
#include <protocols/simple_filters/ResidueCountFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class ResidueCountFilter : public filters::Filter
{
public:
	//default ctor
	ResidueCountFilter();

	ResidueCountFilter(
		ResidueCountFilter const & src
	);

	virtual ~ResidueCountFilter();

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
		moves::DataMap &,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const &
	);

	core::Size
	max_residue_count() const;

	void
	max_residue_count(
		core::Size value
	);

	bool
	enable_max_residue_count() const;

	void
	enable_max_residue_count(
		bool value
	);

	core::Size
	min_residue_count() const;

	void
	min_residue_count(
		core::Size value
	);

	bool
	enable_min_residue_count() const;

	void
	enable_min_residue_count(
		bool value
	);

private:
	core::Size max_residue_count_;
	bool enable_max_residue_count_;
	core::Size min_residue_count_;
	bool enable_min_residue_count_;

};

}
}

#endif
