// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/PdbDataFeatures.hh
/// @author Sam DeLuca
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_PdbDataFeatures_HH
#define INCLUDED_protocols_features_PdbDataFeatures_HH

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/PdbDataFeatures.fwd.hh>

//External

// Project Headers
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class PdbDataFeatures : public protocols::features::FeaturesReporter {
public:
	PdbDataFeatures();

	PdbDataFeatures(PdbDataFeatures const & src);

	virtual ~PdbDataFeatures();

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	void delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

private:

	void load_residue_pdb_identification(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void insert_residue_pdb_identification_rows(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);

	///@brief load the temperature and occupancy information into the
	///PDBInfo object. Backbone atoms are assigned max_bb_temperature
	///and min_bb_occupancy while sidechain atoms are assigned
	///max_sc_temperature and min_sc_occupancy.

	///Note: The information stored in the residue_pdb_confidence table
	///is at the residue level not atom level. Since the temperature and
	///occupancy is usually at the atom level, writing data to the
	///database and reading it back in will result in a loss of
	///information.
	void load_residue_pdb_confidence(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose);

	void insert_residue_pdb_confidence_rows(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	);



};

} // namespace
} // namespace


#endif // include guard
