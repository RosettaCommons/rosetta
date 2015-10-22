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
/// @author Matthew O'Meara (mattjomeara@gmail.com)

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

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>
// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols {
namespace features {

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

void
PairFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_pairs_table_schema(db_session);
}

void
PairFeatures::write_residue_pairs_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column res1_10A_neighbors("res1_10A_neighbors", DbDataTypeOP( new DbInteger() ));
	Column res2_10A_neighbors("res2_10A_neighbors", DbDataTypeOP( new DbInteger() ));
	Column actcoord_dist("actcoord_dist", DbDataTypeOP( new DbReal() ));
	Column polymeric_sequence_dist("polymeric_sequence_dist", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum1);
	primary_key_columns.push_back(resNum2);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(resNum1);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(resNum2);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	GreaterThanConstraintOP res1_10A_neighbors_is_positive( new GreaterThanConstraint(res1_10A_neighbors, 1.0) );
	GreaterThanConstraintOP res2_10A_neighbors_is_positive( new GreaterThanConstraint(res2_10A_neighbors, 1.0) );

	GreaterThanConstraintOP actcoord_dist_is_nonnegative( new GreaterThanConstraint(actcoord_dist, 0) );

	Schema table("residue_pairs", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_constraint(res1_10A_neighbors_is_positive);
	table.add_constraint(res2_10A_neighbors_is_positive);
	table.add_column(res1_10A_neighbors);
	table.add_column(res2_10A_neighbors);
	table.add_column(actcoord_dist);
	table.add_column(polymeric_sequence_dist);

	table.write(db_session);
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
	StructureID const struct_id,
	sessionOP db_session
){
	report_residue_pairs(pose, relevant_residues, struct_id, db_session);
	return 0;
}

void
PairFeatures::report_residue_pairs(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
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

	std::string statement_string = "INSERT INTO residue_pairs (struct_id, resNum1, resNum2, res1_10A_neighbors, res2_10A_neighbors, actcoord_dist, polymeric_sequence_dist) VALUES (?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size resNum1=1; resNum1 <= pose.total_residue(); ++resNum1 ) {
		Residue const & res1( pose.residue(resNum1) );

		Size res1_10A_neighbors(
			tenA.get_node(resNum1)->num_neighbors_counting_self_static());

		// TODO: just iterate over the neighbors of res1
		for ( Size resNum2=resNum1+1; resNum2 <= pose.total_residue(); ++resNum2 ) {
			if ( !check_relevant_residues( relevant_residues, resNum1, resNum2 ) ) continue;
			Residue const & res2( pose.residue(resNum2) );

			Distance const actcoord_dist( res1.actcoord().distance( res2.actcoord() ) );

			if ( actcoord_dist > 10 ) continue;

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
