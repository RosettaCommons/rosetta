// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/StrandHelixGeometryFilter.hh
/// @brief Another filter used in Marcos & Basanta et al. 2017 that needs to be updated.
/// @author Benjamin Basanta (basantab@uw.edu)

#ifndef INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_hh
#define INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_hh

// Unit headers
#include <protocols/fldsgn/filters/StrandHelixGeometryFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh

namespace protocols {
namespace fldsgn {
namespace filters {

///@brief Another filter used in Marcos & Basanta et al. 2017 that needs to be updated.
class StrandHelixGeometryFilter : public protocols::filters::Filter {

public:
	StrandHelixGeometryFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~StrandHelixGeometryFilter() override;

	/// @brief miniimum angle for filtering
	void secstruct( std::string const & ss );

	/// @brief miniimum angle for filtering
	void filter_min_orthoangle( core::Real const r );

	/// @brief maximum angle for filtering
	void filter_max_orthoangle( core::Real const r );

	/// @brief miniimum angle for filtering
	void filter_min_planeangle( core::Real const r );

	/// @brief maximum angle for filtering
	void filter_max_planeangle( core::Real const r );

	/// @brief miniimum angle for filtering
	void filter_min_dist( core::Real const r );

	/// @brief maximum angle for filtering
	void filter_max_dist( core::Real const r );

	/// @brief miniimum angle for filtering
	void strand_id1( core::Size const r );

	/// @brief miniimum angle for filtering
	void strand_id2( core::Size const r );

	/// @brief maximum angle for filtering
	void helix_id( core::Size const r );

	/// @brief set output type
	void output_type( std::string const & s );

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

private:
	// @brief if value is empty, dssp will run for ss definition ( default is emptry )
	mutable std::string secstruct_;

	/// @brief filtered min angle between helix and sheet
	core::Real filter_min_orthoangle_;

	/// @brief filtered max angle between helix and sheet
	core::Real filter_max_orthoangle_;

	/// @brief filtered min angle between helix and sheet
	core::Real filter_min_planeangle_;

	/// @brief filtered max angle between helix and sheet
	core::Real filter_max_planeangle_;

	/// @brief filtered min dist between helix and sheet
	core::Real filter_min_dist_;

	/// @brief filtered max dist between helix and sheet
	core::Real filter_max_dist_;

	/// @brief strand id
	core::Size strand_id1_;

	/// @brief strand id
	core::Size strand_id2_;

	/// @brief helix_id
	core::Size helix_id_;

	/// @brief output_type
	std::string output_type_;

	/// @brief output value, the result of filterring calculation
	mutable core::Real output_value_;

};

} //protocols
} //fldsgn
} //filters

#endif //INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_hh
