// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ProteinBackboneAtomAtomPairFeatures.cc
/// @brief  report atom-atom pair distances between atoms in protein backbones to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ProteinBackboneAtomAtomPairFeatures.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <utility/graph/Graph.hh>
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


// External Headers
#include <cppdb/frontend.h>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ProteinBackboneAtomAtomPairFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::string;
using core::chemical::num_canonical_aas;
using core::chemical::AtomIndices;
using core::pose::Pose;
using core::Size;
using core::Distance;
using core::Vector;
using utility::graph::Graph;
using core::conformation::Residue;
using core::scoring::TenANeighborGraph;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

ProteinBackboneAtomAtomPairFeatures::ProteinBackboneAtomAtomPairFeatures(){}

ProteinBackboneAtomAtomPairFeatures::ProteinBackboneAtomAtomPairFeatures( ProteinBackboneAtomAtomPairFeatures const & ) :
	FeaturesReporter()
{}

ProteinBackboneAtomAtomPairFeatures::~ProteinBackboneAtomAtomPairFeatures()= default;

// XRW TEMP string
// XRW TEMP ProteinBackboneAtomAtomPairFeatures::type_name() const { return "ProteinBackboneAtomAtomPairFeatures"; }

void
ProteinBackboneAtomAtomPairFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_protein_backbone_atom_atom_pairs_table_schema(db_session);
}

void
ProteinBackboneAtomAtomPairFeatures::write_protein_backbone_atom_atom_pairs_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum1("resNum1", DbDataTypeOP( new DbInteger() ));
	Column resNum2("resNum2", DbDataTypeOP( new DbInteger() ));
	Column N_N_dist("N_N_dist", DbDataTypeOP( new DbReal() ));
	Column N_Ca_dist("N_Ca_dist", DbDataTypeOP( new DbReal() ));
	Column N_C_dist("N_C_dist", DbDataTypeOP( new DbReal() ));
	Column N_O_dist("N_O_dist", DbDataTypeOP( new DbReal() ));
	Column N_Ha_dist("N_Ha_dist", DbDataTypeOP( new DbReal() ));
	Column Ca_N_dist("Ca_N_dist", DbDataTypeOP( new DbReal() ));
	Column Ca_Ca_dist("Ca_Ca_dist", DbDataTypeOP( new DbReal() ));
	Column Ca_C_dist("Ca_C_dist", DbDataTypeOP( new DbReal() ));
	Column Ca_O_dist("Ca_O_dist", DbDataTypeOP( new DbReal() ));
	Column Ca_Ha_dist("Ca_Ha_dist", DbDataTypeOP( new DbReal() ));
	Column C_N_dist("C_N_dist", DbDataTypeOP( new DbReal() ));
	Column C_Ca_dist("C_Ca_dist", DbDataTypeOP( new DbReal() ));
	Column C_C_dist("C_C_dist", DbDataTypeOP( new DbReal() ));
	Column C_O_dist("C_O_dist", DbDataTypeOP( new DbReal() ));
	Column C_Ha_dist("C_Ha_dist", DbDataTypeOP( new DbReal() ));
	Column O_N_dist("O_N_dist", DbDataTypeOP( new DbReal() ));
	Column O_Ca_dist("O_Ca_dist", DbDataTypeOP( new DbReal() ));
	Column O_C_dist("O_C_dist", DbDataTypeOP( new DbReal() ));
	Column O_O_dist("O_O_dist", DbDataTypeOP( new DbReal() ));
	Column O_Ha_dist("O_Ha_dist", DbDataTypeOP( new DbReal() ));
	Column Ha_N_dist("Ha_N_dist", DbDataTypeOP( new DbReal() ));
	Column Ha_Ca_dist("Ha_Ca_dist", DbDataTypeOP( new DbReal() ));
	Column Ha_C_dist("Ha_C_dist", DbDataTypeOP( new DbReal() ));
	Column Ha_O_dist("Ha_O_dist", DbDataTypeOP( new DbReal() ));
	Column Ha_Ha_dist("Ha_Ha_dist", DbDataTypeOP( new DbReal() ));

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

	Schema table("protein_backbone_atom_atom_pairs", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);
	table.add_column(N_N_dist);
	table.add_column(N_Ca_dist);
	table.add_column(N_C_dist);
	table.add_column(N_O_dist);
	table.add_column(N_Ha_dist);
	table.add_column(Ca_N_dist);
	table.add_column(Ca_Ca_dist);
	table.add_column(Ca_C_dist);
	table.add_column(Ca_O_dist);
	table.add_column(Ca_Ha_dist);
	table.add_column(C_N_dist);
	table.add_column(C_Ca_dist);
	table.add_column(C_C_dist);
	table.add_column(C_O_dist);
	table.add_column(C_Ha_dist);
	table.add_column(O_N_dist);
	table.add_column(O_Ca_dist);
	table.add_column(O_C_dist);
	table.add_column(O_O_dist);
	table.add_column(O_Ha_dist);
	table.add_column(Ha_N_dist);
	table.add_column(Ha_Ca_dist);
	table.add_column(Ha_C_dist);
	table.add_column(Ha_O_dist);
	table.add_column(Ha_Ha_dist);

	table.write(db_session);
}

utility::vector1<std::string>
ProteinBackboneAtomAtomPairFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

/// @details These atom-atom pairs follow the analysis done in:
///
///  Song Y, Tyka M, Leaver-Fay A, Thompson J, Baker D. Structure
///  guided forcefield optimization. Proteins: Structure, Function,
///  and Bioinformatics. 2011:n/a-n/a. Available at:
///  http://doi.wiley.com/10.1002/prot.23013 [Accessed April 4, 2011].
///
/// The HBond geometries are recoded in the HBondFeatures reporter
Size
ProteinBackboneAtomAtomPairFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	TenANeighborGraph const & tenA(pose.energies().tenA_neighbor_graph());

	std::string statement_string = "INSERT INTO protein_backbone_atom_atom_pairs (struct_id, resNum1, resNum2, N_N_dist, N_Ca_dist, N_C_dist, N_O_dist, N_Ha_dist, Ca_N_dist, Ca_Ca_dist, Ca_C_dist, Ca_O_dist, Ca_Ha_dist, C_N_dist, C_Ca_dist, C_C_dist, C_O_dist, C_Ha_dist, O_N_dist, O_Ca_dist, O_C_dist, O_O_dist, O_Ha_dist, Ha_N_dist, Ha_Ca_dist, Ha_C_dist, Ha_O_dist, Ha_Ha_dist) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size resNum1=1; resNum1 <= pose.size(); ++resNum1 ) {
		Residue const & res1 = pose.residue(resNum1);
		if ( !res1.is_protein() ) continue;

		Vector const & N1(res1.xyz("N"));
		Vector const & Ca1(res1.xyz("CA"));
		Vector const & C1(res1.xyz("C"));
		Vector const & O1(res1.xyz("O"));
		// For glysine, use the hydrogen that is in the same chiral
		// position as HA on other amino acids
		Vector const & HA1(
			res1.aa() != core::chemical::aa_gly ? res1.xyz("HA") : res1.xyz("2HA"));

		for ( Graph::EdgeListConstIter
				ir  = tenA.get_node( resNum1 )->const_edge_list_begin(),
				ire = tenA.get_node( resNum1 )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size resNum2( (*ir)->get_other_ind(resNum1) );

			if ( !check_relevant_residues( relevant_residues, resNum1, resNum2 ) || (resNum1 >= resNum2) ) continue;
			Residue const & res2 = pose.residue(resNum2);
			if ( !res2.is_protein() ) continue;

			Vector const & N2(res2.xyz("N"));
			Vector const & Ca2(res2.xyz("CA"));
			Vector const & C2(res2.xyz("C"));
			Vector const & O2(res2.xyz("O"));
			// For glysine, use the hydrogen that is in the same chiral
			// position as HA on other amino acids
			Vector const & HA2(
				res2.aa() != core::chemical::aa_gly ? res2.xyz("HA") : res2.xyz("2HA"));

			stmt.bind(1,struct_id);
			stmt.bind(2,resNum1);
			stmt.bind(3,resNum2);
			stmt.bind(4,N1.distance(N2));
			stmt.bind(5,N1.distance(Ca2));
			stmt.bind(6,N1.distance(C2));
			stmt.bind(7,N1.distance(O2));
			stmt.bind(8,N1.distance(HA2));
			stmt.bind(9,Ca1.distance(N2));
			stmt.bind(10,Ca1.distance(Ca2));
			stmt.bind(11,Ca1.distance(C2));
			stmt.bind(12,Ca1.distance(O2));
			stmt.bind(13,Ca1.distance(HA2));
			stmt.bind(14,C1.distance(N2));
			stmt.bind(15,C1.distance(Ca2));
			stmt.bind(16,C1.distance(C2));
			stmt.bind(17,C1.distance(O2));
			stmt.bind(18,C1.distance(HA2));
			stmt.bind(19,O1.distance(N2));
			stmt.bind(20,O1.distance(Ca2));
			stmt.bind(21,O1.distance(C2));
			stmt.bind(22,O1.distance(O2));
			stmt.bind(23,O1.distance(HA2));
			stmt.bind(24,HA1.distance(N2));
			stmt.bind(25,HA1.distance(Ca2));
			stmt.bind(26,HA1.distance(C2));
			stmt.bind(27,HA1.distance(O2));
			stmt.bind(28,HA1.distance(HA2));


			basic::database::safely_write_to_database(stmt);
		} //res2
	} //res1
	return 0;
}

std::string ProteinBackboneAtomAtomPairFeatures::type_name() const {
	return class_name();
}

std::string ProteinBackboneAtomAtomPairFeatures::class_name() {
	return "ProteinBackboneAtomAtomPairFeatures";
}

void ProteinBackboneAtomAtomPairFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Records all the distances between backbone atoms -- all possible pairs! -- for each possible residue pair, as features.", attlist );
}

std::string ProteinBackboneAtomAtomPairFeaturesCreator::type_name() const {
	return ProteinBackboneAtomAtomPairFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ProteinBackboneAtomAtomPairFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ProteinBackboneAtomAtomPairFeatures );
}

void ProteinBackboneAtomAtomPairFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinBackboneAtomAtomPairFeatures::provide_xml_schema( xsd );
}

} // namesapce
} // namespace
