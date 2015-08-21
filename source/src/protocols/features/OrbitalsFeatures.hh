// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/OrbitalsFeatures.hh
/// @brief  report Orbital geometry and scores to features Statistics Scientific Benchmark
/// @author Steven Combs

#ifndef INCLUDED_protocols_features_OrbitalsFeatures_hh
#define INCLUDED_protocols_features_OrbitalsFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/OrbitalsFeatures.fwd.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.fwd.hh>

//External

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class OrbitalsFeatures : public protocols::features::FeaturesReporter {
public:
	OrbitalsFeatures();

	OrbitalsFeatures( OrbitalsFeatures const & src );

	virtual ~OrbitalsFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	/// @brief generate the HPOL_orbital table schema
	void
	write_HPOL_orbital_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the HARO_orbital table schema
	void
	write_HARO_orbital_table_schema(
		utility::sql_database::sessionOP db_session) const;

	/// @brief generate the HPOL_orbital table schema
	void
	write_orbital_orbital_table_schema(
		utility::sql_database::sessionOP db_session) const;


public:
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
		utility::sql_database::sessionOP db_session
	);

	void
	report_hpol_orbital_interactions(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session
	);
	void
	report_haro_orbital_interactions(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID const struct_id,
		utility::sql_database::sessionOP db_session
	);
	void
	set_OrbH_features_data(
		core::conformation::Residue  & res1,
		core::conformation::Residue  & res2,
		core::Size const Aindex,
		core::Size const Hindex,
		core::Size const Orbindex,
		numeric::xyzVector<core::Real> const Orbxyz,
		core::Size & resNum2,
		std::string & orbName1,
		std::string & htype2,
		std::string & res2name,
		core::Size & orbNum1,
		core::Size & hpolNum2,
		core::Real & cosAOH,
		core::Real & cosDHO,
		core::Real & chiBDHO,
		core::Real & chiBAOH,
		core::Real & AOH_angle,
		core::Real & DHO_angle,
		core::Real & chiBAHD,
		core::Real & cosAHD,
		core::Real & OrbHdist);
	void
	set_OrbOrb_features_data(
		core::conformation::Residue  & res1,
		core::conformation::Residue  & res2,
		core::Size Aindex,
		core::Size Dindex,
		core::Size Orbindex1,
		core::Size Orbindex2,
		numeric::xyzVector<core::Real> const & Orbxyz1,
		numeric::xyzVector<core::Real> const & Orbxyz2,
		core::Size & resNum2,
		std::string & orbName1,
		std::string & res2name,
		std::string & OrbName2,
		core::Size & orbNum1,
		core::Size & OrbNum2,
		core::Real & cosAOO,
		core::Real & cosDOO,
		core::Real & chiBAOO,
		core::Real & chiBDOO,
		core::Real & AOO_angle,
		core::Real & DOO_angle,
		core::Real & OrbHdist,
		core::Real & cosAOD,
		core::Real & AOD_angle,
		core::Real & chiBAHD,
		core::Real & cosAHD
	);

};

} // features namespace
} // protocols namespace

#endif // include guard
