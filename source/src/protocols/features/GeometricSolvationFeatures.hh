// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/GeometricSolvationFeatures.hh
/// @brief  report hydrogen bonding based solvation model to a features database
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_features_GeometricSolvationFeatures_hh
#define INCLUDED_protocols_features_GeometricSolvationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/GeometricSolvationFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class GeometricSolvationFeatures : public FeaturesReporter {
public:
	GeometricSolvationFeatures(core::scoring::methods::EnergyMethodOptions const& options);

	GeometricSolvationFeatures(GeometricSolvationFeatures const & src);

	~GeometricSolvationFeatures() override= default;

	/// @brief return string with class name
	std::string
	type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const override;

private:
	/// @brief generate the atom_in_residue_pairs table schema
	void
	write_geometric_solvation_table_schema(
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
	core::scoring::geometric_solvation::ExactOccludedHbondSolEnergy geo_sol_energy_;

};


} // features namespace
} // protocols namespace

#endif // include guard
