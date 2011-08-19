// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueSecondaryStructureFeatures.cc
/// @brief  report ResidueSecondaryStructure geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueSecondaryStructureFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <protocols/jumping/Dssp.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using protocols::jumping::Dssp;
using core::Size;
using core::conformation::Residue;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using basic::Tracer;

static Tracer TR("protocols.features.ResidueSecondaryStructureFeatures");

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures() {}

ResidueSecondaryStructureFeatures::ResidueSecondaryStructureFeatures(ResidueSecondaryStructureFeatures const &) :
	FeaturesReporter()
{}

ResidueSecondaryStructureFeatures::~ResidueSecondaryStructureFeatures() {}

string
ResidueSecondaryStructureFeatures::type_name() const { return "ResidueSecondaryStructureFeatures"; }

string
ResidueSecondaryStructureFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS dssp_codes(\n"
		"	code TEXT,\n"
		"	label TEXT,\n"
		"	PRIMARY KEY(code));\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('H', 'H: a-Helix');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('E', 'E: b-Sheet');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('T', 'T: HB Turn');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('G', 'G: 3/10 Helix');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('B', 'B: b-Bridge');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('S', 'S: Bend');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES('I', 'I: pi-Helix');\n"
		"INSERT OR IGNORE INTO dssp_codes VALUES(' ', 'Irregular');\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS residue_secondary_structure(\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	dssp TEXT,\n"
		"	FOREIGN KEY(struct_id, resNum)\n"
		"		REFERENCES residues(struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY(dssp)\n"
		"		REFERENCES dssp_codes(code)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum));\n";
}

Size
ResidueSecondaryStructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	// compute dssp
	protocols::jumping::Dssp all_dssp(pose);

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;

		// don't just push the char into the table.  It will end up an int!
		stringstream dssp; dssp << all_dssp.get_dssp_secstruct(resNum);
		statement stmt = (*db_session)
			<< "INSERT INTO residue_secondary_structure VALUES (?,?,?);"
			<< struct_id
			<< resNum
			<< dssp.str();
		stmt.exec();
	}
	return 0;
}


} // namesapce
} // namespace
