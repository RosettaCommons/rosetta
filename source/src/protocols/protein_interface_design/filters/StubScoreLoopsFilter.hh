// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ContingentFilter.hh
/// @brief A filter that is contingent on some other mover to set its pass/fail value
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_StubScoreLoopsFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_StubScoreLoopsFilter_hh


// Project Headers
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/filters/StubScoreLoopsFilter.fwd.hh>

#include <protocols/constraint_filters/ConstraintScoreCutoffFilter.hh>

#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>


// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class StubScoreLoopsFilter : public protocols::constraint_filters::ConstraintScoreCutoffFilter
{
private:
	typedef protocols::constraint_filters::ConstraintScoreCutoffFilter Parent;
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
public:
	/// @brief default ctor
	StubScoreLoopsFilter();

	~StubScoreLoopsFilter() override;

	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	using protocols::constraint_filters::ConstraintScoreCutoffFilter::get_score; // For the other overloads
	core::Real get_score( core::pose::Pose const & pose ) const override;

private:
	hotspot_hashing::HotspotStubSetOP stub_set_;
	std::string resfile_;
	core::Real cb_force_;
	core::Size loop_start_;
	core::Size loop_stop_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_StubScoreLoopsFilter_HH_
