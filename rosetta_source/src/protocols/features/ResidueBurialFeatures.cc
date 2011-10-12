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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/nv/NVscore.hh>
#include <core/scoring/nv/NVscoreCreator.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

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
		"	struct_id TEXT,\n"
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
Size
ResidueBurialFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
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

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Residue const & res = pose.residue(resNum);

		Size const ten_a_neighbors(tenA.get_node(resNum)->num_neighbors_counting_self_static());
		Size const twelve_a_neighbors(twelveA.get_node(resNum)->num_neighbors_counting_self_static());

		EnergyMap nv_emap;
		nv_score_->residue_energy(res, pose, nv_emap);

		while(true)
		{
			try
			{
				statement stmt = (*db_session)
					<< "INSERT INTO residue_burial VALUES (?,?,?,?,?,?,?,?);"
					<< struct_id
					<< resNum
					<< ten_a_neighbors
					<< twelve_a_neighbors
					<< nv_emap[neigh_vect_raw]
					<< residue_sasa_s[resNum]
					<< residue_sasa_m[resNum]
					<< residue_sasa_l[resNum];
				stmt.exec();
				break;
			}catch(cppdb::cppdb_error &)
			{
				usleep(10);
				continue;
			}
		}
	}
	return 0;
}

} // namespace
} // namespace
