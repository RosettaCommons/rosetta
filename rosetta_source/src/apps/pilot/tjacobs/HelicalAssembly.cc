// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelicalAssembly.cc
///
/// @brief

/// @author Tim Jacobs

// Unit headers
#include <devel/helixAssembly/HelixAssemblyGraph.hh>

// Core headers
#include <core/types.hh>

// Devel headers
#include <devel/init.hh>

//Utility headers
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>


// C++ headers
#include <stdio.h>
#include <time.h>
#include <string>

static basic::Tracer TR("HelicalAssembly");

//double diffclock(clock_t clock1,clock_t clock2)
//{
//	double diffticks=clock1-clock2;
//	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
//	return diffms;
//}

// run protocol

namespace HelicalAssembly {
	basic::options::IntegerOptionKey const path_size( "path_size" );
	basic::options::RealOptionKey const max_rmsd( "max_rmsd" );
	basic::options::RealOptionKey const min_third_helix_rmsd( "min_third_helix_rmsd" );
	basic::options::BooleanOptionKey const simple_print( "simple_print" );
}

int
main( int argc, char * argv [] )
{

	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	option.add( HelicalAssembly::path_size, "");
	option.add( HelicalAssembly::max_rmsd, "Maximum RMSD for a pair of bundles");
	option.add( HelicalAssembly::min_third_helix_rmsd, "");
	option.add( HelicalAssembly::simple_print, "");

	// initialize core
	devel::init(argc, argv);

	core::Real max_rmsd = option[HelicalAssembly::max_rmsd].def(0.7);
	core::Real min_third_helix_rmsd = option[HelicalAssembly::min_third_helix_rmsd].def(5);
	core::Size path_size = option[HelicalAssembly::path_size].def(6);

	if(path_size % 2 != 0){
		utility_exit_with_message("ERROR: Path size must be an even number");
	}

	bool simple_print = option[HelicalAssembly::simple_print].def(false);

	// Initialize DB
	utility::sql_database::sessionOP db_session(basic::database::get_db_session());

	std::string inter_bundle_select_string = "SELECT * FROM pair_comparisons WHERE rmsd <= "+to_string(max_rmsd)+"AND third_helix_rmsd > "+to_string(min_third_helix_rmsd)+";";
	cppdb::statement inter_bundle_select_stmt = basic::database::safely_prepare_statement(inter_bundle_select_string, db_session);
	cppdb::result res1 = basic::database::safely_read_from_database(inter_bundle_select_stmt);

	HelixAssemblyGraph graph(db_session, simple_print);
	while(res1.next()){
		core::Size pair_id_1, pair_id_2;
		core::Real rmsd;
		res1 >> pair_id_1 >> pair_id_2 >> rmsd;

		Node node1(pair_id_1);
		Node node2(pair_id_2);

		graph.addInterBundleEdge(node1, node2);
		graph.addInterBundleEdge(node2, node1);
	}

	std::string intra_bundle_select_string =
	"SELECT p1.pair_id AS pair_id_1, p2.pair_id AS pair_id_2, p3.pair_id AS pair_id_3\n"
	"FROM bundle_pairs p1\n"
	"JOIN bundle_pairs p2 ON\n"
	"	p1.bundle_id=p2.bundle_id AND\n"
	"	p1.pair_id < p2.pair_id\n"
	"JOIN bundle_pairs p3 ON\n"
	"	p2.bundle_id=p3.bundle_id AND\n"
	"	p2.pair_id < p3.pair_id;";

	cppdb::statement intra_bundle_select_stmt = basic::database::safely_prepare_statement(intra_bundle_select_string, db_session);
	cppdb::result res2 = basic::database::safely_read_from_database(intra_bundle_select_stmt);

	while(res2.next()){
		core::Size pair_id_1, pair_id_2, pair_id_3;
		res2 >> pair_id_1 >> pair_id_2 >> pair_id_3;

		Node node1(pair_id_1);
		Node node2(pair_id_2);
		Node node3(pair_id_3);

		graph.addIntraBundleEdge(node1, node2);
		graph.addIntraBundleEdge(node2, node1);

		graph.addIntraBundleEdge(node1, node3);
		graph.addIntraBundleEdge(node3, node1);

		graph.addIntraBundleEdge(node2, node3);
		graph.addIntraBundleEdge(node3, node2);
	}

	//loop through each node to use it as the start node
	std::set<Node> nodes = graph.nodes();

	TR << "Total nodes in graph: " << nodes.size() << std::endl;
	TR << "Total inter-bundle edges: " << graph.total_inter_bundle_edges() << std::endl;
	for(std::set<Node>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
//		clock_t begin=clock();
		utility::vector1<Node> visited;
		visited.push_back((*it));
		graph.depthFirst(visited, path_size);
		clock_t end=clock();
//		cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< endl;
	}
}
