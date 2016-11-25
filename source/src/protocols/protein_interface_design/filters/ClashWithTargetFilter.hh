// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/ClashWithTargetFilter.hh
/// @brief Clash check with Target using motifgrafting
/// @author Lei Shi (shileiustc@gmail.com)

#ifndef INCLUDED_protocols_protein_interface_design_filters_ClashWithTargetFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_ClashWithTargetFilter_hh


#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class ClashWithTargetFilter : public protocols::filters::Filter
{
public:
	ClashWithTargetFilter();
	bool apply( core::pose::Pose const & pose ) const override;
	protocols::filters::FilterOP clone() const override;
	protocols::filters::FilterOP fresh_instance() const override {
		return protocols::filters::FilterOP( new ClashWithTargetFilter() );
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Size compute( core::pose::Pose const & pose ) const;
	virtual ~ClashWithTargetFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & , protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string align_to_pdbname_;
	std::string context_pdbname_;
	core::Size clash_score_cutoff_;
	std::string clash_residues_;
	core::Size ref_start_;
	core::Size ref_end_;
	core::Size pose_start_;
	core::Size pose_end_;
};

} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_ClashWithTargetFilter_HH_

