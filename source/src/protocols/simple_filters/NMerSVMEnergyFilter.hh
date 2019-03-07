// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NMerSVMEnergyFilter.hh
/// @brief definition of filter class NMerSVMEnergyFilter.
/// @author Indigo King (indigo.c.king@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_NMerSVMEnergyFilter_hh
#define INCLUDED_protocols_simple_filters_NMerSVMEnergyFilter_hh

//unit headers
#include <protocols/simple_filters/NMerSVMEnergyFilter.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/NMerSVMEnergy.hh>

namespace protocols {
namespace simple_filters {

class NMerSVMEnergyFilter : public filters::Filter
{
public:
	//default ctor
	NMerSVMEnergyFilter();
	//full ctor
	NMerSVMEnergyFilter(
		core::Real const score_type_threshold,
		std::string string_resnums
	);
	bool apply( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return utility::pointer::make_shared< NMerSVMEnergyFilter >( *this );
	}
	filters::FilterOP fresh_instance() const override{
		return utility::pointer::make_shared< NMerSVMEnergyFilter >();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	// core::Real compute_residue( core::pose::Pose const & pose, core::Size const seqpos, core::Real &, utility::vector1< core::Real > & ) const;
	core::Real compute( core::pose::Pose const &pose ) const;
	~NMerSVMEnergyFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real score_type_threshold_;
	bool dump_table_;
	std::string string_resnums_;
	core::scoring::methods::EnergyMethodOptions opts_;
	core::scoring::methods::NMerSVMEnergy energy_method_ = core::scoring::methods::NMerSVMEnergy( opts_ );
	bool count_eps_;
	bool rank_as_score_;
	core::Real ep_cutoff_;
};

}
}

#endif
