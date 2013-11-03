// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/TotalSasaFilter.hh
/// @brief definition of filter class TotalSasaFilter.
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_TotalSasaFilter_hh
#define INCLUDED_protocols_simple_filters_TotalSasaFilter_hh

#include <protocols/simple_filters/TotalSasaFilter.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_filters {

class TotalSasaFilter : public filters::Filter
{
public:
	TotalSasaFilter();
	TotalSasaFilter( core::Real const lower_threshold, bool const hydrophobic=false, bool const polar=false, core::Real const upper_threshold=100000000.0 );

	bool apply( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	core::pack::task::TaskFactoryOP task_factory();
	/// @brief Set the task factory to limit SASA calculation to packable residues
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	filters::FilterOP clone() const;
	filters::FilterOP fresh_instance() const;

	virtual ~TotalSasaFilter();

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	core::Real lower_threshold_;
	core::Real upper_threshold_;
	core::pack::task::TaskFactoryOP taskfactory_; /// Residue subset definition. Only the surface area from packable residues.
	bool hydrophobic_, polar_; /// count only hydrophobics? polars?
};

}
}

#endif
