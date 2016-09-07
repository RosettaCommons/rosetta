// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/ResidueGridScoresFeatures.hh
/// @brief detailed per atom scores of Scoring Grids
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_ResidueGridScoresFeatures_hh
#define INCLUDED_protocols_features_ResidueGridScoresFeatures_hh

#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueGridScoresFeatures.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>


namespace protocols {
namespace features {

class ResidueGridScoresFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueGridScoresFeatures();

	ResidueGridScoresFeatures( ResidueGridScoresFeatures const & src );

	~ResidueGridScoresFeatures() override;

	/// @brief return string with class name
	std::string
	type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/) override;

private:
	char chain_;


};

}
}


#endif
