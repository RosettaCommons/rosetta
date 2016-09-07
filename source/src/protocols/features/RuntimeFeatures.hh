// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_features_RuntimeFeatures_HH
#define INCLUDED_protocols_features_RuntimeFeatures_HH

// Unit headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/RuntimeFeatures.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.hh>

// RosettaScripts headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace features {

class RuntimeFeatures : public FeaturesReporter {

public:

	/// @brief Default constructor.
	RuntimeFeatures();

	/// @brief Default destructor.
	~RuntimeFeatures() override;

	/// @copydoc FeaturesReporter::get_name
	std::string type_name() const override { return "RuntimeFeatures"; }

	/// @brief The runtime features links to a structure ID, so the
	/// StructureFeatures reporter must also be present.
	utility::vector1<std::string> features_reporter_dependencies() const override;

	/// @copydoc FeaturesReporter::write_schema_to_db
	void write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

	/// @brief Report runtime information for the current job.
	core::Size report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session) override;

};

}
}

#endif
