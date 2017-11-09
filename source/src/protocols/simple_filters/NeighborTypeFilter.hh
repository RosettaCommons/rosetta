// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NeighborTypeFilter.hh
/// @brief definition of filter class InterfaceSasaFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_NeighborTypeFilter_hh
#define INCLUDED_protocols_simple_filters_NeighborTypeFilter_hh

#include <protocols/simple_filters/NeighborTypeFilter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

class NeighborTypeFilter : public filters::Filter
{
public:
	NeighborTypeFilter() : filters::Filter( "NeighborType" ){};
	NeighborTypeFilter( std::string const & target_residue,
		utility::vector1< bool > const & residue_types,
		core::Real const distance_threshold
	) :
		filters::Filter( "NeighborType" ),
		target_residue_(target_residue),
		residue_types_(residue_types),
		distance_threshold_(distance_threshold)
	{}

	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	std::vector< core::Size > compute( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const override {
		return filters::FilterOP( new NeighborTypeFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override{
		return filters::FilterOP( new NeighborTypeFilter() );
	}
	void clear() override { residue_types_.clear(); }
	~NeighborTypeFilter() override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string target_residue_;
	utility::vector1< bool > residue_types_;
	core::Real distance_threshold_;
};

}
}

#endif
