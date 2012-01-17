// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/PairFeatures.cc
/// @brief  report orbital geometry and scores to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/PairFeatures.hh>

// Project Headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::chemical::AtomIndices;
using core::pose::Pose;
using core::Size;
using core::Distance;
using core::Vector;
using core::graph::Graph;
using core::conformation::Residue;
using core::scoring::TenANeighborGraph;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;

PairFeatures::PairFeatures(){}

PairFeatures::PairFeatures( PairFeatures const & ) :
	FeaturesReporter()
{}

PairFeatures::~PairFeatures(){}

string
PairFeatures::type_name() const { return "PairFeatures"; }

string
PairFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS residue_pairs (\n"
		"	struct_id TEXT,\n"
		"	resNum1 INTEGER,\n"
		"	resNum2 INTEGER,\n"
		"	res1_10A_neighbors INTEGER,\n"
		"	res2_10A_neighbors INTEGER,\n"
		"	actcoord_dist REAL,\n"
		"	polymeric_sequence_dist INTEGER,\n"
		"	FOREIGN KEY (struct_id, resNum1)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	FOREIGN KEY (struct_id, resNum2)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	CONSTRAINT res1_10A_neighbors_is_positive CHECK (res1_10A_neighbors >= 1),\n"
		"	CONSTRAINT res2_10A_neighbors_is_positive CHECK (res2_10A_neighbors >= 1),\n"
		"	CONSTRAINT actcoord_dist_is_nonnegative CHECK (actcoord_dist >= 0),\n"
		"	PRIMARY KEY(struct_id, resNum1, resNum2));";
}

utility::vector1<std::string>
PairFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}


Size
PairFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){
	report_residue_pairs(pose, relevant_residues, struct_id, db_session);
	return 0;
}

void
PairFeatures::report_residue_pairs(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	// assert pose.update_residue_neighbors() has been called:
	runtime_assert(
		 !pose.conformation().structure_moved() &&
		 pose.energies().residue_neighbors_updated());

	// I would like to assert that the actcoords are up to date but that
	// isn't currently possible.
	//pose.update_actcoords();

	TenANeighborGraph const & tenA( pose.energies().tenA_neighbor_graph() );

	std::string statement_string = "INSERT INTO residue_pairs VALUES (?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size resNum1=1; resNum1 <= pose.total_residue(); ++resNum1){
		if(!relevant_residues[resNum1]) continue;
		Residue res1( pose.residue(resNum1) );

		Size res1_10A_neighbors(
			tenA.get_node(resNum1)->num_neighbors_counting_self_static());

		// TODO: just iterate over the neighbors of res1
		for(Size resNum2=resNum1+1; resNum2 <= pose.total_residue(); ++resNum2){
			if(!relevant_residues[resNum2]) continue;
			Residue res2( pose.residue(resNum2) );

			Distance const actcoord_dist( res1.actcoord().distance( res2.actcoord() ) );

			if( actcoord_dist > 10 ) continue;

			Size res2_10A_neighbors(
				tenA.get_node(resNum2)->num_neighbors_counting_self_static());

			int polymeric_sequence_dist(res1.polymeric_sequence_distance(res2));

			stmt.bind(1,struct_id);
			stmt.bind(2,resNum1);
			stmt.bind(3,resNum2);
			stmt.bind(4,res1_10A_neighbors);
			stmt.bind(5,res2_10A_neighbors);
			stmt.bind(6,actcoord_dist);
			stmt.bind(7,polymeric_sequence_dist);
			basic::database::safely_write_to_database(stmt);
		}
	}
}

} // namesapce
} // namespace
