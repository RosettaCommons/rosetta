// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/DdGFeatures.hh
/// @brief  report the per-residue ddG score to the features database
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_DdGFeatures_hh
#define INCLUDED_protocols_features_DdGFeatures_hh

// Unit Headers
#include <protocols/features/DdGFeatures.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

#include <protocols/simple_filters/DdGScan.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace features {

class DdGFeatures : public protocols::features::FeaturesReporter {
public:
	DdGFeatures();

	DdGFeatures(
		protocols::simple_filters::DdGScanOP ddG_scan_mover
	);

	DdGFeatures(DdGFeatures const & src);

	virtual ~DdGFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db (
		utility::sql_database::sessionOP db_session
	) const;

	void
	write_ddG_table_schema (
		utility::sql_database::sessionOP db_session
	) const;

private:
	/// @brief generate the residue_total_scores_1b table schema
	void
	insert_ddG_rows (
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session
	) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	void
	parse_my_tag (
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
	);


	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

private:

	protocols::simple_filters::DdGScanOP ddG_scan_mover_;

};

} // namespace
} // namespace

#endif // include guard
