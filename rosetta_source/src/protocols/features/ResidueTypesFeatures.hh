// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueTypesFeatures.hh
/// @brief  report ResidueTypes features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ResidueTypesFeatures_hh
#define INCLUDED_protocols_features_ResidueTypesFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/ResidueTypesFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <string>

#include <core/chemical/ResidueType.fwd.hh>


namespace protocols{
namespace features{

class ResidueTypesFeatures : public protocols::features::FeaturesReporter {
public:
	ResidueTypesFeatures();

	/// Undefined, commenting out to fix PyRosetta build  ResidueTypesFeatures(core::scoring::ScoreFunctionOP scfxn);

	/// Undefined, commenting out to fix PyRosetta build  ResidueTypesFeatures( ResidueTypesFeatures const & src );

	virtual ~ResidueTypesFeatures();

	///@breif return string with class name
	std::string
	type_name() const;

	///@breif return sql statements that setup the right tables
	std::string
	schema() const;

	///@breif collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

private:

	void
	report_residue_type(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_atom(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_bond(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_cut_bond(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_chi(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_chi_rotamer(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_proton_chi(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_properties(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

	void
	report_residue_type_variant(
		std::string const & residue_type_set_name,
		core::chemical::ResidueType const & res_type,
		utility::sql_database::sessionOP db_session) const;

private:

	// should match version string in residue type parameter sets
	core::Real const version_;

};

} // namespace
} // namespace

#endif // include guard
