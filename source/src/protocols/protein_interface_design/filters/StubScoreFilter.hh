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

#ifndef INCLUDED_protocols_protein_interface_design_filters_StubScoreFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_StubScoreFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilter.fwd.hh>

#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <utility/vector1.hh>


// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace filters {

class StubScoreFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
public:
	/// @brief default ctor
	StubScoreFilter();
	/// @brief Constructor with a single target residue
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~StubScoreFilter();
	void stub_sets( utility::vector1< StubSetStubPos > const & stub_sets );
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	utility::vector1< StubSetStubPos > stub_sets_;
	core::Size host_chain_;
	core::Real cb_force_;
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_StubScoreFilter_HH_
