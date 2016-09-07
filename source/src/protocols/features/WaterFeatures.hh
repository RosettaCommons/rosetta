// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/WaterFeatures.hh
/// @brief  report geometry between water and hbond sites
/// @author Kevin Houlihan

#ifndef INCLUDED_protocols_features_WaterFeatures_hh
#define INCLUDED_protocols_features_WaterFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/WaterFeatures.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class WaterFeatures : public FeaturesReporter {
public:
	WaterFeatures();

	WaterFeatures(core::Real acc_dist_cutoff, core::Real don_dist_cutoff);

	WaterFeatures(WaterFeatures const & src);

	~WaterFeatures() override= default;

	core::Real
	acc_dist_cutoff() const { return acc_dist_cutoff_; }

	core::Real
	don_dist_cutoff() const { return don_dist_cutoff_; }

	void
	acc_dist_cutoff(core::Real d) { acc_dist_cutoff_ = d; }

	void
	don_dist_cutoff(core::Real d) { don_dist_cutoff_ = d; }

	/// @brief return string with class name
	std::string
	type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

private:
	/// @brief generate the salt_bridges table schema
	void
	write_water_hbond_geom_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const &,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

private:

	core::Real acc_dist_cutoff_;
	core::Real don_dist_cutoff_;
	core::Real ahd_cutoff_;

	utility::vector1< std::pair<std::string, std::string> > names_for_water_;
};


} // features namespace
} // protocols namespace

#endif // include guard
