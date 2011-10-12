// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinBackboneTorsionAngleFeatures.cc
/// @brief  report Backbone Torsional Angle features
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinBackboneTorsionAngleFeatures.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// Platform Headers
#include <core/pose/Pose.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using cppdb::statement;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;


ProteinBackboneTorsionAngleFeatures::ProteinBackboneTorsionAngleFeatures(){}

ProteinBackboneTorsionAngleFeatures::ProteinBackboneTorsionAngleFeatures( ProteinBackboneTorsionAngleFeatures const & ) :
	FeaturesReporter()
{}

ProteinBackboneTorsionAngleFeatures::~ProteinBackboneTorsionAngleFeatures(){}

string
ProteinBackboneTorsionAngleFeatures::type_name() const { return "ProteinBackboneTorsionAngleFeatures"; }

string
ProteinBackboneTorsionAngleFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS protein_backbone_torsion_angles (\n"
		"	struct_id TEXT,\n"
		"	resNum INTEGER,\n"
		"	phi REAL,\n"
		"	psi REAL,\n"
		"	omega REAL,\n"
		"	FOREIGN KEY (struct_id, resNum)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY (struct_id, resNum));";
}

Size
ProteinBackboneTorsionAngleFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!relevant_residues[i]) continue;

		Residue const & resi = pose.residue(i);
		Real phi  (resi.mainchain_torsion(1));
		Real psi  (resi.mainchain_torsion(2));
		Real omega(resi.mainchain_torsion(3));
		while(true)
		{
			try
			{
				statement stmt = (*db_session) <<
					"INSERT INTO protein_backbone_torsion_angles VALUES (?,?,?,?,?)" <<
					struct_id << i << phi << psi << omega;
				stmt.exec();
				break;
			}catch(cppdb::cppdb_error &)
			{
				usleep(10);
				continue;
			}
		}
	}
	return 0;
}

} // namesapce
} // namespace
