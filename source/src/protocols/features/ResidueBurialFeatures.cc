// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/ResidueBurialFeatures.cc
/// @brief  report residue burial to features statistics scientific benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/features/ResidueBurialFeatures.hh>

// Project Headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/nv/NVscore.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/ResidueBurialFeaturesCreator.hh>


namespace protocols {
namespace features {

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using core::id::AtomID_Map;
using core::scoring::calc_per_atom_sasa;
using core::scoring::EnergyMap;
using core::scoring::TenANeighborGraph;
using core::scoring::TwelveANeighborGraph;
using core::scoring::neigh_vect_raw;
using core::scoring::methods::EnergyMethodOptions;
using core::scoring::nv::NVscore;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

ResidueBurialFeatures::ResidueBurialFeatures() :
	nv_score_(core::scoring::nv::NVscoreOP( new NVscore() ))
{}

ResidueBurialFeatures::ResidueBurialFeatures(ResidueBurialFeatures const & src) :
	FeaturesReporter(),
	nv_score_(src.nv_score_)
{}

ResidueBurialFeatures::~ResidueBurialFeatures()= default;

// XRW TEMP string
// XRW TEMP ResidueBurialFeatures::type_name() const { return "ResidueBurialFeatures"; }

void
ResidueBurialFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_residue_burial_table_schema(db_session);
}

void
ResidueBurialFeatures::write_residue_burial_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column resNum("resNum", DbDataTypeOP( new DbInteger() ));
	Column ten_a_neighbors("ten_a_neighbors", DbDataTypeOP( new DbInteger() ));
	Column twelve_a_neighbors("twelve_a_neighbors", DbDataTypeOP( new DbInteger() ));
	Column neigh_vect_raw("neigh_vect_raw", DbDataTypeOP( new DbReal() ));
	Column sasa_r100("sasa_r100", DbDataTypeOP( new DbReal() ));
	Column sasa_r140("sasa_r140", DbDataTypeOP( new DbReal() ));
	Column sasa_r200("sasa_r200", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(resNum);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(resNum);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);

	Schema table("residue_burial", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(ten_a_neighbors);
	table.add_column(twelve_a_neighbors);
	table.add_column(neigh_vect_raw);
	table.add_column(sasa_r100);
	table.add_column(sasa_r140);
	table.add_column(sasa_r200);

	table.write(db_session);
}

utility::vector1<std::string>
ResidueBurialFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

Size
ResidueBurialFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){

	TenANeighborGraph const & tenA( pose.energies().tenA_neighbor_graph() );
	TwelveANeighborGraph const & twelveA( pose.energies().twelveA_neighbor_graph() );

	Real const probe_radius_s(1.0);
	AtomID_Map< Real > atom_sasa_s;
	vector1< Real > residue_sasa_s;
	calc_per_atom_sasa( pose, atom_sasa_s, residue_sasa_s, probe_radius_s);

	Real const probe_radius_m(1.4);
	AtomID_Map< Real > atom_sasa_m;
	vector1< Real > residue_sasa_m;
	calc_per_atom_sasa( pose, atom_sasa_m, residue_sasa_m, probe_radius_m);

	Real const probe_radius_l(2.0);
	AtomID_Map< Real > atom_sasa_l;
	vector1< Real > residue_sasa_l;
	calc_per_atom_sasa( pose, atom_sasa_l, residue_sasa_l, probe_radius_l);

	std::string statement_string = "INSERT INTO residue_burial (struct_id, resNum, ten_a_neighbors, twelve_a_neighbors, neigh_vect_raw, sasa_r100, sasa_r140, sasa_r200) VALUES (?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for ( Size resNum=1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues( relevant_residues, resNum ) ) continue;
		Residue const & res = pose.residue(resNum);

		Size const ten_a_neighbors(tenA.get_node(resNum)->num_neighbors_counting_self_static());
		Size const twelve_a_neighbors(twelveA.get_node(resNum)->num_neighbors_counting_self_static());

		EnergyMap nv_emap;
		nv_score_->residue_energy(res, pose, nv_emap);

		stmt.bind(1,struct_id);
		stmt.bind(2,resNum);
		stmt.bind(3,ten_a_neighbors);
		stmt.bind(4,twelve_a_neighbors);
		stmt.bind(5,nv_emap[neigh_vect_raw]);
		stmt.bind(6,residue_sasa_s[resNum]);
		stmt.bind(7,residue_sasa_m[resNum]);
		stmt.bind(8,residue_sasa_l[resNum]);
		basic::database::safely_write_to_database(stmt);
	}
	return 0;
}

std::string ResidueBurialFeatures::type_name() const {
	return class_name();
}

std::string ResidueBurialFeatures::class_name() {
	return "ResidueBurialFeatures";
}

void ResidueBurialFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Record the number of neighbors at different radii and sasa with different probe radii", attlist );
}

std::string ResidueBurialFeaturesCreator::type_name() const {
	return ResidueBurialFeatures::class_name();
}

protocols::features::FeaturesReporterOP
ResidueBurialFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new ResidueBurialFeatures );
}

void ResidueBurialFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueBurialFeatures::provide_xml_schema( xsd );
}


} // namespace
} // namespace
