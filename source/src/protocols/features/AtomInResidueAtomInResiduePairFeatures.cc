// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/AtomInResidueAtomInResiduePairFeatures.cc
/// @brief  report atom-atom pair geometry and scores to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/AtomInResidueAtomInResiduePairFeatures.hh>

// Project Headers
#include <core/chemical/AA.hh>
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
#include <basic/database/schema_generator/DbDataType.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray5D.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

//Auto Headers
namespace protocols {
namespace features {

using std::string;
using core::chemical::num_canonical_aas;
using core::chemical::AtomIndices;
using core::pose::Pose;
using core::Size;
using core::Distance;
using core::Vector;
using core::graph::Graph;
using core::conformation::Residue;
using core::scoring::TenANeighborGraph;
using ObjexxFCL::FArray5D;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

AtomInResidueAtomInResiduePairFeatures::AtomInResidueAtomInResiduePairFeatures(){}

AtomInResidueAtomInResiduePairFeatures::AtomInResidueAtomInResiduePairFeatures( AtomInResidueAtomInResiduePairFeatures const & ) :
	FeaturesReporter()
{}

AtomInResidueAtomInResiduePairFeatures::~AtomInResidueAtomInResiduePairFeatures()= default;

string
AtomInResidueAtomInResiduePairFeatures::type_name() const { return "AtomInResidueAtomInResiduePairFeatures"; }

void
AtomInResidueAtomInResiduePairFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_atom_in_residue_pairs_table_schema(db_session);
}

void
AtomInResidueAtomInResiduePairFeatures::write_atom_in_residue_pairs_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column residue_type1("residue_type1", DbDataTypeOP( new DbText() ));
	Column atom_type1("atom_type1", DbDataTypeOP( new DbText() ));
	Column residue_type2("residue_type2", DbDataTypeOP( new DbText() ));
	Column atom_type2("atom_type2", DbDataTypeOP( new DbText() ));
	Column distance_bin("distance_bin", DbDataTypeOP( new DbText() ));
	Column count("count", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_type1);
	primary_key_columns.push_back(atom_type1);
	primary_key_columns.push_back(residue_type2);
	primary_key_columns.push_back(atom_type2);
	primary_key_columns.push_back(distance_bin);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);

	GreaterThanConstraintOP count_is_non_negative( new GreaterThanConstraint(count, 0) );

	Schema table("atom_in_residue_pairs", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_constraint(count_is_non_negative);
	table.add_column(count);

	table.write(db_session);
}


utility::vector1<std::string>
AtomInResidueAtomInResiduePairFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
AtomInResidueAtomInResiduePairFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	report_atom_pairs(pose, relevant_residues, struct_id, db_session);
	return 0;
}

/// @detail This is very similar in spirit to the potential described in
///
///﻿Lu H, Skolnick J. A distance-dependent atomic knowledge-based potential for improved protein structure selection. Proteins. 2001;44(3):223-32. Available at: http://www.ncbi.nlm.nih.gov/pubmed/11455595.
///
/// However, they use different distance bins.  Here, [0,1), ...,
/// [9,10) are used because they are easy and as they report the the
/// paper, most of the signal comes in the 3.5-6.5 range.  To get the
/// molar fraction of atom types--since the types are unique within
/// each residue type, there is exactly one per residue of that type.
/// Therefore this information can be extracted from the Residues
/// table when needed.  It may make sense to include it here if it
/// turns to to be too cumbersom to get those quantities.
///
/// TODO: Expand for not just canonical residue types


void
AtomInResidueAtomInResiduePairFeatures::report_atom_pairs(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){

	// assert pose.update_residue_neighbors() has been called:
	runtime_assert(
		!pose.conformation().structure_moved() &&
		pose.energies().residue_neighbors_updated());

	Size const max_res(num_canonical_aas);
	Size const max_atm(30); // check this
	Size const dist_bins(15);
	FArray5D< Size > counts;
	counts.dimension(max_res, max_atm, max_res, max_atm, dist_bins, 1);

	TenANeighborGraph const & tenA( pose.energies().tenA_neighbor_graph() );


	for ( Size resNum1=1; resNum1 <= pose.total_residue(); ++resNum1 ) {
		Residue res1( pose.residue(resNum1) );

		for ( Graph::EdgeListConstIter
				ir  = tenA.get_node( resNum1 )->const_edge_list_begin(),
				ire = tenA.get_node( resNum1 )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size resNum2( (*ir)->get_other_ind(resNum1) );
			if ( !check_relevant_residues( relevant_residues, resNum1, resNum2 ) ) continue;

			Residue res2( pose.residue(resNum2) );

			for ( Size atmNum1=1; atmNum1 <= res1.natoms(); ++atmNum1 ) {
				Vector const & atm1_xyz( res1.xyz(atmNum1) );

				for ( Size atmNum2=1; atmNum2 <= res2.natoms(); ++atmNum2 ) {
					Vector const & atm2_xyz( res2.xyz(atmNum2) );

					Size const dist_bin(static_cast<Size>(ceil(atm1_xyz.distance(atm2_xyz))));
					if ( dist_bin < 15 ) {
						counts(res1.aa(), atmNum1, res2.aa(), atmNum2, dist_bin) += 1;
					}
				}
			}
		}
	}

	std::string stmt_string = "INSERT INTO atom_in_residue_pairs (struct_id, residue_type1, atom_type1, residue_type2, atom_type2, distance_bin, count) VALUES (?,?,?,?,?,?,?);";
	cppdb::statement stmt(basic::database::safely_prepare_statement(stmt_string,db_session));

	for ( Size aa1=1; aa1 <= max_res; ++aa1 ) {
		for ( Size aa2=1; aa2 <= max_res; ++aa2 ) {
			for ( Size atmNum1=1; atmNum1 <= max_atm; ++atmNum1 ) {
				for ( Size atmNum2=1; atmNum2 <= max_atm; ++atmNum2 ) {
					for ( Size dist_bin=1; dist_bin <= 15; ++dist_bin ) {
						Size const count(counts(aa1, atmNum1, aa2, atmNum2, dist_bin));
						stmt.bind(1,struct_id);
						stmt.bind(2,aa1);
						stmt.bind(3,atmNum1);
						stmt.bind(4,aa2);
						stmt.bind(5,atmNum2);
						stmt.bind(6,dist_bin);
						stmt.bind(7,count);
						basic::database::safely_write_to_database(stmt);
					}
				}
			}
		}
	}
}

} // namesapce
} // namespace
