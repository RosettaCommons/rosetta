// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ResidueBurialFeatures.cc
/// @brief  report residue burial to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ResidueBurialFeatures.hh>

// Project Headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/nv/NVscore.hh>
// AUTO-REMOVED #include <core/scoring/nv/NVscoreCreator.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>


namespace protocols{
namespace features{

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using core::id::AtomID_Map;
using core::scoring::calc_per_atom_sasa;
using core::scoring::EnergyMap;
using core::scoring::ScoreFunctionOP;
using core::scoring::TenANeighborGraph;
using core::scoring::TwelveANeighborGraph;
using core::scoring::neigh_vect_raw;
using core::scoring::methods::EnergyMethodOptions;
using core::scoring::nv::NVscore;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

ResidueBurialFeatures::ResidueBurialFeatures() :
	nv_score_(new NVscore())
{}

ResidueBurialFeatures::ResidueBurialFeatures(ResidueBurialFeatures const & src) :
	FeaturesReporter(),
	nv_score_(src.nv_score_)
{}

ResidueBurialFeatures::~ResidueBurialFeatures(){}

string
ResidueBurialFeatures::type_name() const { return "ResidueBurialFeatures"; }

string
ResidueBurialFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS residue_burial (\n"
		"	struct_id BLOB,\n"
		"	resNum INTEGER,\n"
		"	ten_a_neighbors INTEGER,\n"
		"	twelve_a_neighbors INTEGER,\n"
		"	neigh_vect_raw REAL,\n"
		"	sasa_r100 REAL,\n"
		"	sasa_r140 REAL,\n"
		"	sasa_r200 REAL,\n"
		"	FOREIGN KEY (struct_id, resNum)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY (struct_id, resNum));\n";
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
	boost::uuids::uuid const struct_id,
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

	std::string statement_string = "INSERT INTO residue_burial VALUES (?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
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

} // namespace
} // namespace
