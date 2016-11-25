// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/filters/DisulfideFilter.hh
/// @brief Filters for interfaces which could form a disulfide bond between
/// docking partners.
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_AtomicContactFilter_hh
#define INCLUDED_protocols_simple_filters_AtomicContactFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/ResId.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

/// @brief detects atomic (<4Ang) contacts between any two atoms of two residues
class AtomicContactFilter : public protocols::filters::Filter, public protocols::moves::ResId
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	AtomicContactFilter();
	AtomicContactFilter( core::Size const res1, core::Size const res2, core::Real const distance=4.0, bool const sidechain=1, bool const backbone=0, bool const protons=0 );
	bool apply( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new AtomicContactFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new AtomicContactFilter() );
	}

	~AtomicContactFilter() override= default;
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
	core::Size residue1_;/*, residue2_ residue2 is managed by the ResId baseclass*/
	utility::vector1< core::Size > range1_;
	core::Real distance_;
	bool sidechain_, backbone_, protons_;
};

} // filters
} // protocols

#endif
