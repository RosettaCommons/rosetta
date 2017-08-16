// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueTypesFeatures.hh
/// @brief  report ResidueTypes features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_features_ResidueTypesFeatures_hh
#define INCLUDED_protocols_features_ResidueTypesFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueTypesFeatures.fwd.hh>
#include <core/chemical/ResidueDatabaseIO.hh>

//External

// Project Headers
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <string>

#include <core/chemical/ResidueType.fwd.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace features {

class ResidueTypesFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueTypesFeatures();

	/// Undefined, commenting out to fix PyRosetta build  ResidueTypesFeatures(core::scoring::ScoreFunctionOP scfxn);

	/// Undefined, commenting out to fix PyRosetta build  ResidueTypesFeatures( ResidueTypesFeatures const & src );

	~ResidueTypesFeatures() override;

	/// @brief return string with class name

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

	std::string
	type_name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	// should match version string in residue type parameter sets
	core::Real version_;
	core::chemical::ResidueDatabaseIO residue_dbio_;

};

} // namespace
} // namespace

#endif // include guard
