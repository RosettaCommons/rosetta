// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinBackboneAtomAtomPairFeatures.cc
/// @brief  report atom-atom pair distances between atoms in protein backbones to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinBackboneAtomAtomPairFeatures.hh>

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

// ObjexxFCL Headers
#include <ObjexxFCL/FArray5D.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

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
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

ProteinBackboneAtomAtomPairFeatures::ProteinBackboneAtomAtomPairFeatures(){}

ProteinBackboneAtomAtomPairFeatures::ProteinBackboneAtomAtomPairFeatures( ProteinBackboneAtomAtomPairFeatures const & ) :
	FeaturesReporter()
{}

ProteinBackboneAtomAtomPairFeatures::~ProteinBackboneAtomAtomPairFeatures(){}

string
ProteinBackboneAtomAtomPairFeatures::type_name() const { return "ProteinBackboneAtomAtomPairFeatures"; }

string
ProteinBackboneAtomAtomPairFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS protein_backbone_atom_atom_pairs (\n"
		"	struct_id INTEGER,\n"
		"	resNum1 TEXT,\n"
		"	resNum2 TEXT,\n"
		"	N_N_dist REAL,\n"
		"	N_Ca_dist REAL,\n"
		"	N_C_dist REAL,\n"
		"	N_O_dist REAL,\n"
		"	Ca_N_dist REAL,\n"
		"	Ca_Ca_dist REAL,\n"
		"	Ca_C_dist REAL,\n"
		"	Ca_O_dist REAL,\n"
		"	C_N_dist REAL,\n"
		"	C_Ca_dist REAL,\n"
		"	C_C_dist REAL,\n"
		"	C_O_dist REAL,\n"
		"	O_N_dist REAL,\n"
		"	O_Ca_dist REAL,\n"
		"	O_C_dist REAL,\n"
		"	O_O_dist REAL,\n"
		"	FOREIGN KEY (struct_id)\n"
		"		REFERENCES structures (struct_id)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY (struct_id, resNum1, resNum2));\n";
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
	Size const struct_id,
	sessionOP db_session
){
	TenANeighborGraph const & tenA(pose.energies().tenA_neighbor_graph());

	for(Size resNum1=1; resNum1 <= pose.total_residue(); ++resNum1){
		if(!relevant_residues[resNum1]) continue;
		Residue const & res1 = pose.residue(resNum1);
		if(!res1.is_protein()) continue;

		Vector const & N1(res1.xyz("N"));
		Vector const & Ca1(res1.xyz("CA"));
		Vector const & C1(res1.xyz("C"));
		Vector const & O1(res1.xyz("O"));

		for ( Graph::EdgeListConstIter
			ir  = tenA.get_node( resNum1 )->const_edge_list_begin(),
			ire = tenA.get_node( resNum1 )->const_edge_list_end();
			ir != ire; ++ir ) {
			Size resNum2( (*ir)->get_other_ind(resNum1) );
			if(!relevant_residues[resNum2] || (resNum1 >= resNum2)) continue;
			Residue const & res2 = pose.residue(resNum2);
			if(!res2.is_protein()) continue;

			Vector const & N2(res2.xyz("N"));
			Vector const & Ca2(res2.xyz("CA"));
			Vector const & C2(res2.xyz("C"));
			Vector const & O2(res2.xyz("O"));

			statement stmt = (*db_session)
				<< "INSERT INTO protein_backbone_atom_atom_pairs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
				<< struct_id << resNum1 << resNum2
				<< N1.distance(N2)
				<< N1.distance(Ca2)
				<< N1.distance(C2)
				<< N1.distance(O2)
				<< Ca1.distance(N2)
				<< Ca1.distance(Ca2)
				<< Ca1.distance(C2)
				<< Ca1.distance(O2)
				<< C1.distance(N2)
				<< C1.distance(Ca2)
				<< C1.distance(C2)
				<< C1.distance(O2)
				<< O1.distance(N2)
				<< O1.distance(Ca2)
				<< O1.distance(C2)
				<< O1.distance(O2);
			stmt.exec();
		} //res2
	} //res1
	return 0;
}
} // namesapce
} // namespace
