// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_features_TotalScoreFeatures_HH
#define INCLUDED_protocols_features_TotalScoreFeatures_HH

// Unit headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/TotalScoreFeatures.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace features {

class TotalScoreFeatures : public FeaturesReporter {

public:

	/// @brief Default constructor.
	TotalScoreFeatures();

	/// @brief Constructor with score function argument.
	TotalScoreFeatures(core::scoring::ScoreFunctionOP scorefxn);

	/// @brief Default destructor.
	~TotalScoreFeatures();

	/// @copydoc FeaturesReporter::type_name
	std::string type_name() const;

	/// @brief Get the score function being reported by this object.
	core::scoring::ScoreFunctionCOP scorefxn() const;

	/// @brief Set the score function being reported by this object.
	void scorefxn(core::scoring::ScoreFunctionOP scorefxn);

	/// @copydoc FeaturesReporter::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose);

	/// @copydoc FeaturesReporter::write_schema_to_db
	void write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	/// @copydoc FeaturesReporter::report_features
	core::Size report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

private:
	core::scoring::ScoreFunctionOP scorefxn_;

};

}
}

#endif
