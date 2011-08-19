// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinRMSDFeatures.cc
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinRMSDFeatures.hh>

// Project Headers
#include <core/pose/util.hh>

// Platform Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>


// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <string>
#include <map>
#include <list>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::list;
using core::Size;
using core::pose::Pose;
using core::pose::PoseCOP;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;

string
ProteinRMSDFeatures::type_name() const { return "ProteinRMSDFeatures"; }


ProteinRMSDFeatures::ProteinRMSDFeatures(
	PoseCOP reference_pose ) :
	reference_pose_(reference_pose)
{}

string
ProteinRMSDFeatures::schema() const {
	string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3") {
		return
			"CREATE TABLE IF NOT EXISTS protein_rmsd (\n"
			"	struct_id INTEGER,\n"
			"	reference_tag TEXT,\n"
			"	protein_CA REAL,\n"
			"	protein_CA_or_CB REAL,\n"
			"	protein_backbone REAL,\n"
			"	protein_backbone_including_O REAL,\n"
			"	protein_backbone_sidechain_heavyatom REAL,\n"
			"	heavyatom REAL,\n"
			"	nbr_atom REAL,\n"
			"	all_atom REAL,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, reference_tag));";
	}else if(db_mode == "mysql") {
		return
			"CREATE TABLE IF NOT EXISTS protein_rmsd (\n"
			"	struct_id INTEGER,\n"
			"	reference_tag VARCHAR(255),\n"
			"	protein_CA DOUBLE,\n"
			"	protein_CA_or_CB DOUBLE,\n"
			"	protein_backbone DOUBLE,\n"
			"	protein_backbone_including_O DOUBLE,\n"
			"	protein_backbone_sidechain_heavyatom DOUBLE,\n"
			"	heavyatom DOUBLE,\n"
			"	nbr_atom DOUBLE,\n"
			"	all_atom DOUBLE,\n"
			"	FOREIGN KEY (struct_id)\n"
			"		REFERENCES structures (struct_id),\n"
			"	PRIMARY KEY(struct_id, reference_tag));";
	}else {
		return "";
	}
}


Size
ProteinRMSDFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	using namespace core::scoring;

	list< Size > subset_residues;
	for(Size i = 1; i <= relevant_residues.size(); ++i){
		if(relevant_residues[i]) subset_residues.push_back(i);
	}


	statement stmt = (*db_session)
		<< "INSERT INTO protein_rmsd VALUES (?,?,?,?,?,?,?,?,?,?);"
		<< struct_id
		<< find_tag(*reference_pose_)
		<< rmsd_with_super(*reference_pose_, pose, subset_residues, is_protein_CA)
		<< rmsd_with_super(
			*reference_pose_, pose, subset_residues, is_protein_CA_or_CB)
		<< rmsd_with_super(
			*reference_pose_, pose, subset_residues, is_protein_backbone)
		<< rmsd_with_super(
			*reference_pose_, pose, subset_residues, is_protein_backbone_including_O)
		<< rmsd_with_super(
			*reference_pose_, pose, subset_residues, is_protein_sidechain_heavyatom)
		<< rmsd_with_super(*reference_pose_, pose, subset_residues, is_heavyatom)
		<< rmsd_with_super(*reference_pose_, pose, subset_residues, is_nbr_atom)
		<< all_atom_rmsd(*reference_pose_, pose, subset_residues);
	stmt.exec();


	return 0;
}

} //namesapce
} //namespace
