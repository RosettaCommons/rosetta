// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/evolution/AASynthesisFitnessCostFilter.hh

#ifndef INCLUDED_protocols_evolution_AASynthesisFitnessCostFilter_hh
#define INCLUDED_protocols_evolution_AASynthesisFitnessCostFilter_hh

//unit headers
#include <protocols/evolution/AASynthesisFitnessCostFilter.fwd.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace evolution {

/// @brief test whether a pose contains a comment that evaluates to a predefined value. This is useful in controlling execution flow in RosettaScripts.
class AASynthesisFitnessCost : public filters::Filter
{
public:
	AASynthesisFitnessCost();
	~AASynthesisFitnessCost() override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AASynthesisFitnessCost( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new AASynthesisFitnessCost() );
	}

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const &) override;
	core::Real compute( core::pose::Pose const & pose ) const;

	core::Size threshold() const { return threshold_; }
	void threshold( core::Size const & t ) { threshold_ = t; }

	core::Real fitnessCostPerATP() const { return fitnessCostPerATP_; }
	void fitnessCostPerATP( core::Real const & c ) { fitnessCostPerATP_ = c; }

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Size threshold_; //dflt ""; define the comment value
	core::Real fitnessCostPerATP_;
};
}
}

#endif
