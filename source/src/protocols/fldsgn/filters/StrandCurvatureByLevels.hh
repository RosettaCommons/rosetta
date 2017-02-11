// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/StrandCurvatureByLevels.hh
/// @brief Newer version of filter used in Marcos & Basanta et al. 2017
/// @author Benjamin Basanta (basantab@uw.edu)

#ifndef INCLUDED_protocols_fldsgn_filters_StrandCurvatureByLevels_hh
#define INCLUDED_protocols_fldsgn_filters_StrandCurvatureByLevels_hh

// Unit headers
#include <protocols/fldsgn/filters/StrandCurvatureByLevels.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>

// Utility headers

// Parser headers

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh
#include <utility/vector1.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>


namespace protocols {
namespace fldsgn {
namespace filters {

///@brief Newer version of filter used in Marcos & Basanta et al. 2017
class StrandCurvatureByLevels : public protocols::filters::Filter {

public:

	typedef core::Real Real;
	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef std::string String;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;

	StrandCurvatureByLevels();

	// destructor (important for properly forward-declaring smart-pointer members)
	~StrandCurvatureByLevels() override;

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

	// @brief set filtered sheet_topology by SrandPairingSetOP
	//void filtered_sheet_topology( String const & sheet_topology );

	/// @brief miniimum angle for filtering
	void secstruct( String const & ss );

	void bend_level( Size const r );

	/// @brief miniimum angle for filtering
	void filter_min_bend( Real const r );

	/// @brief maximum angle for filtering
	void filter_max_bend( Real const r );

	void twist_level( Size const r );

	/// @brief miniimum twist for filtering
	void filter_min_twist( Real const r );

	/// @brief maximum twist for filtering
	void filter_max_twist( Real const r );

	/// @brief miniimum angle for filtering
	void strand_id( Size const r );

	/// @brief set output type
	void output_type( String const & s );

	void concavity_reference_residue( String const & ss );

	void concavity_direction( bool const & r );

public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	virtual  core::Real compute( Pose const & pose ) const;

private:

	/// @brief if value is empty, dssp will run for ss definition ( default is emptry )
	mutable String secstruct_;

	Size bend_level_;

	/// @brief filtered min angle between helix and sheet
	Real filter_min_bend_;

	/// @brief filtered max angle between helix and sheet
	Real filter_max_bend_;

	Size twist_level_;

	/// @brief filtered min twist between helix and sheet
	Real filter_min_twist_;
	/// @brief filtered max twist between helix and sheet
	Real filter_max_twist_;

	/// @brief strand id
	Size strand_id_;

	/// @brief output_type
	String output_type_;

	/// @brief output value, the result of filterring calculation
	mutable Real output_value_;

	String concavity_reference_residue_;

	bool concavity_direction_;
};

} //protocols
} //fldsgn
} //filters

#endif //INCLUDED_protocols_fldsgn_filters_StrandCurvatureByLevels_hh
