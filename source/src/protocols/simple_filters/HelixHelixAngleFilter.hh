// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/HelixHelixAngleFilter.hh
/// @brief  HelixHelixAngleFilter computes wither the crossing angle or distance between two TMs
/// @details HelixHelixAngleFilter computes 1 of 3 parameters: 1. crossing angle 2. atomic distance 3. vector distance
/// 1. the crossing angle at the closest point between two vectors representing the two helices
/// 2. the shortest distance between two atoms on either TM
/// 3. the shortest distance (point of approach) between the vectors representing the two helices
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_simple_filters_HelixHelixAngleFilter_hh
#define INCLUDED_protocols_simple_filters_HelixHelixAngleFilter_hh

#include <protocols/simple_filters/HelixHelixAngleFilter.fwd.hh>


// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/exit.hh>

#include <string>
#include <utility/vector1.hh>

namespace protocols {
namespace simple_filters {

class HelixHelixAngleFilter : public filters::Filter
{
public:
	HelixHelixAngleFilter() : filters::Filter( "HelixHelixAngle" ) {}
	bool apply( core::pose::Pose const & pose ) const override;
	void report( std::ostream & out, core::pose::Pose const & pose ) const override;
	core::Real report_sm( core::pose::Pose const & pose ) const override;
	core::Real compute( core::pose::Pose const & pose ) const;

	// @brief iterate all atoms in both TMs, return the minimal distance between them
	core::Real calc_shortest_dist_by_atoms( core::pose::Pose const & pose  ) const;

	// @brief find point of approach (closest points) on the vectors, and return the distance
	std::pair< numeric::xyzVector< core::Real >, numeric::xyzVector< core::Real > >
	find_closest_pnts( utility::vector1< numeric::xyzVector<core::Real>> l1, utility::vector1< numeric::xyzVector<core::Real>> l2 ) const;

	// @brief claculate a vector representing the TM (a 3D linear regression)
	utility::vector1< numeric::xyzVector < core::Real > >
	find_helix_vector( core::pose::Pose const & pose, core::Size start, core::Size end  ) const;

	filters::FilterOP clone() const override {
		return filters::FilterOP( new HelixHelixAngleFilter( *this ) );
	}
	filters::FilterOP fresh_instance() const override {
		return filters::FilterOP( new HelixHelixAngleFilter() );
	}

	virtual ~HelixHelixAngleFilter();
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose ) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Size start_helix_1_ = 0;
	core::Size start_helix_2_ = 0;
	core::Size end_helix_1_ = 0;
	core::Size end_helix_2_ = 0;
	core::Real angle_min_ = 40.0;
	core::Real angle_max_ = 100.0;
	core::Size helix_num_1_ = 1;
	core::Size helix_num_2_ = 2;
	core::Real dist_min_ = 0.0;
	core::Real dist_max_ = 5.0;
	std::string angle_or_dist_ = "angle";
	bool by_helices_ = false;
	bool dist_by_atom_ = true;

	std::pair< core::Size, core::Size >
	find_closest_res( core::pose::Pose const & pose, core::Size s1, core::Size e1, core::Size s2, core::Size e2 ) const;
};

}
}
#endif
