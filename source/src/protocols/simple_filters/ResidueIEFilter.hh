// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueIEFilter.hh
/// @brief definition of filter class ResidueIEFilter.
/// @author Sagar Khare (khares@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_ResidueIEFilter_hh
#define INCLUDED_protocols_simple_filters_ResidueIEFilter_hh

#include <protocols/simple_filters/ResidueIEFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>
#include <set>

namespace protocols {
namespace simple_filters {

class ResidueIEFilter : public filters::Filter
{
public:
	ResidueIEFilter();

	ResidueIEFilter(
		std::string const & resnums,
		std::string const & restype,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::scoring::ScoreType const score_type = core::scoring::total_score,
		core::Real const threshold = 0.0,
		bool const whole_pose = false,
		bool const whole_interface = false,
		core::Size const rb_jump = 1,
		core::Real const interface_distance_cutoff =  8.0,
		core::Real max_penalty = 1000.0,
		core::Real penalty_factor = 1.0,
		bool const report_energy = false );

	ResidueIEFilter( ResidueIEFilter const &init );
	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new ResidueIEFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new ResidueIEFilter() );
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	~ResidueIEFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	std::string const & resnums() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	core::scoring::ScoreType score_type() const;
	core::Real threshold() const;
	void resnums( std::string const & resnums_str );
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void score_type( core::scoring::ScoreType score_type );
	void threshold( core::Real const th );

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::set< core::Size > compute_resnums( core::pose::Pose const & pose ) const;
	core::Real penalty_from_score( core::Real const res_intE ) const;

private:
	std::string resnum_str_;
	std::string restype_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_pose_;
	bool whole_interface_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;
	core::Real max_penalty_;
	core::Real penalty_factor_;
	bool report_energy_;
	bool use_resE_;
	core::select::residue_selector::ResidueSelectorCOP selector_;

};

}
}

#endif
