// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueFeatures.cc
/// @brief  report residue features to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueFeatures.hh>

// Package Headers
#include <core/scoring/hbonds/HBondSet.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/types.hh>
#include <protocols/jumping/Dssp.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::pose::PoseOP;
using core::conformation::Residue;
using core::scoring::EnergyMap;
using core::scoring::ScoreFunctionOP;
using core::scoring::getScoreFunction;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreTypes;
using numeric::xyzVector;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;


ResidueFeatures::ResidueFeatures() :
	scfxn_(getScoreFunction())
{}

ResidueFeatures::ResidueFeatures(
	ScoreFunctionOP scfxn) :
	scfxn_(scfxn)
{}

ResidueFeatures::ResidueFeatures( ResidueFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ResidueFeatures::~ResidueFeatures()
{}

string
ResidueFeatures::type_name() const { return "ResidueFeatures"; }

string
ResidueFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS residues (\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	name3 TEXT,\n"
		"	res_type TEXT,\n"
		"	FOREIGN KEY (struct_id)\n"
		"		REFERENCES structures (struct_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	CONSTRAINT resNum_is_positive CHECK (resNum >= 1),\n"
		"	PRIMARY KEY(struct_id, resNum));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS residue_scores_1b (\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	score_type TEXT,\n"
		"	score_value REAL,\n"
		"	context_dependent  INTEGER,\n"
		"	FOREIGN KEY (struct_id, resNum)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum, score_type));\n"
		"\n"
		"CREATE TABLE IF NOT EXISTS residue_scores_2b (\n"
		"	struct_id INTEGER,\n"
		"	resNum1 INTEGER,\n"
		"	resNum2 INTEGER,\n"
		"	score_type TEXT,\n"
		"	score_value REAL,\n"
		"	context_dependent  INTEGER,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY(struct_id, resNum1, resNum2, score_type));\n";
}

Size
ResidueFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	insert_residue_rows(pose, relevant_residues, struct_id, db_session);
	insert_residue_scores_rows(pose, relevant_residues, struct_id, db_session );

	return 0;
}


void
ResidueFeatures::insert_residue_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue res = pose.residue(resNum);

		string const name3( res.name3() );
		string const res_type( res.name() );

		statement stmt = (*db_session)
			<< "INSERT INTO residues VALUES (?,?,?,?);"
			<< struct_id
			<< resNum
			<< name3
			<< res_type;
		stmt.exec();
	}
}

///@details
///
/// * Score types are either one body, two body, or whole structure and
/// can either dependend on the context or not.
///
/// * Whole structure terms are reported in the structure_scores table
/// along with the totals from the remaining terms.
///
/// * The one and two body terms are broken up into two different
/// tables because they are parametrized differently (one residue vs
/// two residues).
///
/// * Residues are identified by Rosetta's residue numbering scheme.
/// To convert to the PDB residue numbering scheme, join with the
/// residues_pdb table.
///
/// * Although two body terms can be with in the same residue, these
/// 'intrares' score terms are reported with the two body terms where
/// resNum == otherResNum.
///
/// * Two body terms always are reported so that resNum1 <= resNum2.
///
/// * Values for score terms are only reported if they are non-zero.
///
/// * Values of two body energies are only reported when both residues
/// are true in the relevant_residues vector.

void
ResidueFeatures::insert_residue_scores_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	// I would like to assert that this has been called, but I don't know how
	//scfxn_.setup_for_scoring( pose );

	ScoreTypes ci_1b( scfxn_->ci_1b_types() );
	ScoreTypes cd_1b( scfxn_->cd_1b_types() );
	ScoreTypes ci_2b( scfxn_->ci_2b_types() );
	ScoreTypes cd_2b( scfxn_->cd_2b_types() );


	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue rsd( pose.residue(resNum) );
		{ // Context Independent One Body Energies
			EnergyMap emap;
			scfxn_->eval_ci_1b(rsd, pose, emap);
			for(ScoreTypes::const_iterator st = ci_1b.begin(), ste = ci_1b.end(); st != ste; ++st){
				if(!emap[*st]) continue;

				string const score_type( ScoreTypeManager::name_from_score_type(*st) );
				Real const score_value( emap[*st] );

				statement stmt = (*db_session)
					<< "INSERT INTO residue_scores_1b VALUES (?,?,?,?,?);"
					<< struct_id
					<< resNum
					<< score_type
					<< score_value
					<< false;
				stmt.exec();
			}
		}
		{ // Context Dependent One Body Energies
			EnergyMap emap;
			scfxn_->eval_cd_1b(rsd, pose, emap);
			for(ScoreTypes::const_iterator
				st = cd_1b.begin(), ste = cd_1b.end();
				st != ste; ++st){

				if(!emap[*st]) continue;

				string const score_type( ScoreTypeManager::name_from_score_type(*st) );
				Real const score_value( emap[*st] );
				statement stmt = (*db_session)
					<< "INSERT INTO residue_scores_1b VALUES (?,?,?,?,?);"
					<< struct_id
					<< resNum
					<< score_type
					<< score_value
					<< true;
				stmt.exec();
			}
		}

		// Two Body Energies
		for(Size otherResNum=resNum+1; otherResNum <= pose.total_residue(); ++otherResNum){
			if(!relevant_residues[otherResNum]) continue;
			if(!scfxn_->are_they_neighbors(pose, resNum, otherResNum)) continue;
			Residue otherRsd( pose.residue(otherResNum) );
			{ // Context Independent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_ci_2b(rsd, otherRsd, pose, emap);
				for(ScoreTypes::const_iterator st = ci_2b.begin(), ste = ci_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;

					string const score_type(ScoreTypeManager::name_from_score_type(*st));
					Real const score_value( emap[*st] );

					statement stmt = (*db_session)
						<< "INSERT INTO residue_scores_2b VALUES (?,?,?,?,?,?);"
						<< struct_id
						<< resNum
						<< otherResNum
						<< score_type
						<< score_value
						<< false;
					stmt.exec();
				}
			}
			{ // Context Dependent Two Body Energies
				EnergyMap emap;
				scfxn_->eval_cd_2b(rsd, otherRsd, pose, emap);
				for(ScoreTypes::const_iterator st = cd_2b.begin(), ste = cd_2b.end(); st != ste; ++st){
					if(!emap[*st]) continue;
					string const score_type(ScoreTypeManager::name_from_score_type(*st));
					Real const score_value( emap[*st] );
					statement stmt = (*db_session)
						<< "INSERT INTO residue_scores_2b VALUES (?,?,?,?,?,?);"
						<< struct_id
						<< resNum
						<< otherResNum
						<< score_type
						<< score_value
						<< true;
					stmt.exec();
				}
			}
		} // End Two Body Energies
	} // End res1 for loop
} // End function body


} // namesapce
} // namespace
