// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueBurialFilter.hh
/// @brief definition of filter class ResidueBurialFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_ResidueBurialFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueBurialFilter_hh

#include <protocols/simple_filters/ResidueBurialFilter.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace simple_filters {

class ResidueBurialFilter : public filters::Filter
{
public:
	ResidueBurialFilter();
	// SJF 3Sep13 ResidueBurialFilter( core::Size const target_residue, core::Size const neighbors, core::Real const distance_threshold );
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override;
	filters::FilterOP fresh_instance() const override;
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );

	~ResidueBurialFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	core::Real residue_fraction_buried() const { return residue_fraction_buried_; }
	void residue_fraction_buried( core::Real const r ){ residue_fraction_buried_ = r; }

	std::string residue() const { return residue_; }
	void residue( std::string const & s ){ residue_ = s; }

	core::Size neighbors() const{ return neighbors_; }
	void neighbors( core::Size const n ){ neighbors_ = n; }

	core::Real distance_threshold() const{ return distance_threshold_; }
	void distance_threshold( core::Real const r ){ distance_threshold_ = r; }

private:
	std::string residue_; // dflt ""; save as string, parse the actual residue number at apply
	//SJF 3Sep13 core::Size target_residue_;
	core::Size neighbors_;
	core::Real distance_threshold_;
	core::pack::task::TaskFactoryOP task_factory_; /// used to determine which residues to check for burial dynamically. All designable residues will be checked, and if any of them is buried, returns true
	core::Real residue_fraction_buried_; // dflt 0.0001; what fraction of the residues specified by the task_factory should be buried for the filter to pass

};

}
}

#endif
