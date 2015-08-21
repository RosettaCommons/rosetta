// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueSecondaryStructureFeatures.hh
/// @brief  report ResidueSecondaryStructure geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ResidueSecondaryStructureFeatures_hh
#define INCLUDED_protocols_features_ResidueSecondaryStructureFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueSecondaryStructureFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class ResidueSecondaryStructureFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueSecondaryStructureFeatures();

	ResidueSecondaryStructureFeatures( ResidueSecondaryStructureFeatures const & src );

	virtual ~ResidueSecondaryStructureFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

};

} // namespace
} // namespace

#endif // include guard
