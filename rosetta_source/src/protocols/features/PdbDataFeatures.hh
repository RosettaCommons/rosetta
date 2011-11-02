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

#ifndef INCLUDED_protocols_features_PdbDataFeatures_HH_
#define INCLUDED_protocols_features_PdbDataFeatures_HH_

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/PdbDataFeatures.fwd.hh>

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

	///@brief retrun sql statementes that setup the right tables
	std::string
	schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);

	void delete_record(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

private:

	void insert_residue_pdb_identification_rows(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose);

	void load_residue_pdb_identification(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose);

	void insert_residue_pdb_confidence_rows(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session,
		core::pose::Pose const & pose);

	// don't load residue_pdb_confidence because we're not extracting it
	// at the atomic level

};

}
}


#endif /* PDBDATAFEATURES_HH_ */
