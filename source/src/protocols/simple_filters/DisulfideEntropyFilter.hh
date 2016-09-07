// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/DisulfideEntropyFilter.hh
/// @brief Filter on the entropic effect of disulfide linkage
/// @author Gabriel Rocklin (grocklin@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_DisulfideEntropyFilter_hh
#define INCLUDED_protocols_simple_filters_DisulfideEntropyFilter_hh

//unit headers
#include <protocols/simple_filters/DisulfideEntropyFilter.fwd.hh>

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

class DisulfideEntropyFilter : public filters::Filter
{
public:
	//default ctor
	DisulfideEntropyFilter();

	DisulfideEntropyFilter(core::Real tightness_, core::Real lower_bound_)

	;

	DisulfideEntropyFilter(
		DisulfideEntropyFilter const & src
	);

	~DisulfideEntropyFilter() override;

	bool
	apply(
		core::pose::Pose const & pose
	) const override;

	filters::FilterOP
	clone() const override;

	filters::FilterOP
	fresh_instance() const override;

	void
	report(
		std::ostream & out,
		core::pose::Pose const & pose
	) const override;

	core::Real
	report_sm(
		core::pose::Pose const & pose
	) const override;

	core::Real
	compute_residual(
		core::pose::Pose const & pose
	) const;

	core::Real
	compute(
		core::pose::Pose const & pose
	) const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	core::Real
	lower_bound() const;

	void
	lower_bound(
		core::Real value
	);

	core::Real
	tightness() const;

	void
	tightness(
		core::Real value
	);


private:
	core::Real tightness_;
	core::Real lower_bound_;


};

}
}

#endif
