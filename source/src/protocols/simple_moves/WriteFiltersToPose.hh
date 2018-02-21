// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/WriteFiltersToPose.hh
/// @brief Writes all filter results to the pose file.
/// @author Cody Krivacic (krivacic@berkeley.edu)

#ifndef INCLUDED_protocols_simple_filters_WriteFiltersToPose_hh
#define INCLUDED_protocols_simple_filters_WriteFiltersToPose_hh

//unit headers
#include <protocols/simple_moves/WriteFiltersToPose.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_moves {

/// @brief Writes the results from all filters to the pose file.
/// @details This mover goes through all filters in the RosettaScripts file and writes their
/// user-defined name, (optionally) their filter type, and their numerical result to the
/// output .pdb file. The include_type option writes the filter type.
/// The user can also specify a prefix to append before the user-defined
/// filter name; for instance, let's say there is a PackStat filter, as well as several movers like
/// fastdesign and docking. The "prefix" option allows the user to have several versions of
/// WriteFiltersToPose which take the filter results from different steps in the process.
/// If the RosettaScripts XML file has these instances of WriteFiltersToPose:
///      <WriteFiltersToPose name="writer1" prefix="post_fastdesign_"/>
///      <WriteFiltersToPose name="writer2" prefix="post_docking_" include_type="true"/>
/// If the only filter is a PackStat filter named "packing",
/// the output pdb will have two lines at the end that look something like this:
///      post_docking_packing PackStat 0.6794
///      post_fastdesign_packing 0.6751
/// If there are multiple filters, each instance of this mover will write all
/// filters to the pose.
///
/// Note: If you would like to only write a single filter result to your pose,
/// use FilterReportAsPoseExtraScoresMover

class WriteFiltersToPose : public moves::Mover
{
private:
	protocols::filters::Filters_map filters_;
	std::string prefix_;
	bool include_type_ = false;
public:
	//typedef utility::tag::TagCOP TagCOP;
	//default ctor

	//protocols::filters::Filters_map f;
	WriteFiltersToPose();

	void apply( core::pose::Pose & ) override;

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new WriteFiltersToPose( *this ));
	}

	protocols::moves::MoverOP fresh_instance() const override{
		return protocols::moves::MoverOP( new WriteFiltersToPose() );
	}

	void parse_my_tag(utility::tag::TagCOP,basic::datacache::DataMap &,protocols::filters::Filters_map const &f,protocols::moves::Movers_map const &,core::pose::Pose const &) override;
	std::string get_name() const override;

	static
	std::string
	mover_name();

	static
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};
}
}

#endif
