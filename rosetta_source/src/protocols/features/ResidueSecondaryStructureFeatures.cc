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
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <basic/database/sql_utils.hh>
// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::scoring::dssp::Dssp;
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
		"	PRIMARY KEY(struct_id, resNum));\n"
    
        "CREATE TABLE IF NOT EXISTS helix_segments(\n"
        "	struct_id INTEGER,\n"    
        "	helix_id INTEGER,\n"
        "	residue_begin INTEGER,\n"
        "	residue_end INTEGER,\n"
        "	FOREIGN KEY(struct_id, residue_begin, residue_begin)\n"
        "		REFERENCES residues(struct_id, resNum, resNum)\n"
        "		DEFERRABLE INITIALLY DEFERRED,\n"
        "	PRIMARY KEY(struct_id, helix_id));\n"
    
        "CREATE TABLE IF NOT EXISTS beta_segments(\n"
        "	struct_id INTEGER,\n"    
        "	beta_id INTEGER,\n"
        "	residue_begin INTEGER,\n"
        "	residue_end INTEGER,\n"
        "	FOREIGN KEY(struct_id, residue_begin, residue_begin)\n"
        "		REFERENCES residues(struct_id, resNum, resNum)\n"
        "		DEFERRABLE INITIALLY DEFERRED,\n"
        "	PRIMARY KEY(struct_id, beta_id));\n"
    
//        "CREATE TABLE IF NOT EXISTS loop_segments(\n"
//        "	struct_id INTEGER,\n"    
//        "	loop_id INTEGER,\n"
//        "	residue_begin INTEGER,\n"
//        "	residue_end INTEGER,\n"
//        "	FOREIGN KEY(struct_id, residue_begin, residue_begin)\n"
//        "		REFERENCES residues(struct_id, resNum, resNum)\n"
//        "		DEFERRABLE INITIALLY DEFERRED,\n"
//        "	PRIMARY KEY(struct_id, loop_id));\n"
    ;
}

utility::vector1<std::string>
ResidueSecondaryStructureFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueSecondaryStructureFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	// compute dssp
	core::scoring::dssp::Dssp all_dssp(pose);

    //stores the secondary structure for the current stretch of secondary structure
    string segment_secondary;
    Size segment_begin;
    Size segment_end;
    Size helix_counter=1;
    Size beta_counter=1;
    Size loop_counter=1;
	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;


		if(!pose.residue(resNum).is_protein()){
			// Due to limitations with the current DSSP code,
			// the indexing gets off after a non-protein residue
			// and leads to a segmentation fault.
			break;
		}

        string residue_secondary = string(1, all_dssp.get_dssp_secstruct(resNum));
        if(residue_secondary != segment_secondary){
            
            if(resNum > 1){
                segment_end=resNum-1;
                if(segment_secondary == "H"){
                    std::string statement_string = "INSERT INTO helix_segments VALUES (?,?,?,?);";
                    statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
                    stmt.bind(1,struct_id);
                    stmt.bind(2,helix_counter);
                    stmt.bind(3,segment_begin);
                    stmt.bind(4,segment_end);
                    basic::database::safely_write_to_database(stmt);
                    
                    ++helix_counter;
                }
                else if(segment_secondary == "E"){
                    std::string statement_string = "INSERT INTO beta_segments VALUES (?,?,?,?);";
                    statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
                    stmt.bind(1,struct_id);
                    stmt.bind(2,beta_counter);
                    stmt.bind(3,segment_begin);
                    stmt.bind(4,segment_end);
                    basic::database::safely_write_to_database(stmt);
                    
                    ++beta_counter;
                }
//                else if(segment_secondary == "L"){
//                    std::string statement_string = "INSERT INTO loop_segments VALUES (?,?,?,?);";
//                    statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
//                    stmt.bind(1,struct_id);
//                    stmt.bind(2,loop_counter);
//                    stmt.bind(3,segment_begin);
//                    stmt.bind(4,segment_end);
//                    basic::database::safely_write_to_database(stmt);
//                    
//                    ++loop_counter;
//                }
            }
            
            segment_secondary = residue_secondary;
            segment_begin=resNum;
        }
        
		std::string statement_string = "INSERT INTO residue_secondary_structure VALUES (?,?,?);";
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,residue_secondary);
		basic::database::safely_write_to_database(stmt);

	}
	return 0;
}


} // namesapce
} // namespace
