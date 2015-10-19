// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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
#include <devel/sewing/ProteinAssemblyGraph.hh>
#include <devel/sewing/ThreeHelixGraph.hh>
#include <devel/sewing/ThreeHelixRepeatGraph.hh>
#include <devel/sewing/FourHelixGraph.hh>
#include <devel/sewing/FourHelixRepeatGraph.hh>

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
#include <basic/options/keys/sewing.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <stdio.h>
#include <time.h>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "HelicalAssembly" );

//double diffclock(clock_t clock1,clock_t clock2)
//{
//	double diffticks=clock1-clock2;
//	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
//	return diffms;
//}

// run protocol

namespace HelicalAssembly {
	basic::options::IntegerOptionKey const num_helices( "num_helices" );
	basic::options::RealOptionKey const max_rmsd( "max_rmsd" );
	basic::options::RealOptionKey const min_clash_score( "min_clash_score" );
	basic::options::BooleanOptionKey const simple_print( "simple_print" );
	basic::options::BooleanOptionKey const find_cycles( "find_cycles" );
	basic::options::IntegerOptionKey const node_size( "node_size" );
	basic::options::IntegerOptionKey const repeat_prints( "repeat_prints" );
}

int
main( int argc, char * argv [] )
{
	try {

	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace devel::sewing;

	option.add( HelicalAssembly::num_helices, "Number of helices in final structure");
	option.add( HelicalAssembly::max_rmsd, "Maximum RMSD for a pair of bundles");
	option.add( HelicalAssembly::min_clash_score, "");
	option.add( HelicalAssembly::simple_print, "");
	option.add( HelicalAssembly::find_cycles, "Look for only paths that cycle");
	option.add( HelicalAssembly::node_size, "3 or 4 helix bundles nodes");
	option.add( HelicalAssembly::repeat_prints, "number of repeats to print");

	// initialize core
	devel::init(argc, argv);

	// Initialize DB
	utility::sql_database::sessionOP db_session(
		basic::database::get_db_session());

	core::Real max_rmsd = option[HelicalAssembly::max_rmsd].def(0.7);
	core::Real min_clash_score = option[HelicalAssembly::min_clash_score].def(5);
	core::Size num_helices = option[HelicalAssembly::num_helices].def(4);
	core::Size node_size = option[HelicalAssembly::node_size].def(3);
	if(node_size != 3 && node_size !=4)
	{
		utility_exit_with_message("You must select either 3 or 4 for node size");
	}

	bool simple_print = option[HelicalAssembly::simple_print].def(false);
	bool find_cycles = option[HelicalAssembly::find_cycles].def(false);
	core::Size repeat_prints = option[HelicalAssembly::repeat_prints].def(2);

	ProteinAssemblyGraphOP graph;
	if(find_cycles){
		if(node_size==3){
			graph = new ThreeHelixRepeatGraph(db_session);
		}
		else if(node_size==4){
			graph = new FourHelixRepeatGraph(db_session);
		}
		else{
			utility_exit_with_message("Sewing protocol only currently supports 3 and 4 helix nodes.");
		}
	}
	else{
		if(node_size==3){
			graph = new ThreeHelixGraph(db_session);
		}
		else if(node_size==4){
			graph = new FourHelixGraph(db_session);
		}
		else{
			utility_exit_with_message("Sewing protocol only currently supports 3 and 4 helix nodes.");
		}
	}
	graph->populate_graph(max_rmsd, min_clash_score);

	/****FIND PATHS/TREES****/
	utility::vector1<ProteinAssemblyGraph::EdgeList> trees = graph->find_all_trees(num_helices);
	TR << "Found " << trees.size() << " total trees" << std::endl;

	//remove all duplicates (paths which contain the same set of helices)
	utility::vector1<ProteinAssemblyGraph::EdgeList> trees_no_dups;
	std::set< std::set<core::Size> > helix_id_sets;
	for(core::Size i=1; i<=trees.size(); ++i){

		//create a set of helix ides from the current list of edges (std::pairs)
		std::set<core::Size> helix_ids;
		for(core::Size j=1; j<=trees[i].size(); ++j){
			helix_ids.insert(trees[i][j].first->helix_1_id());
			helix_ids.insert(trees[i][j].first->helix_2_id());

			helix_ids.insert(trees[i][j].second->helix_1_id());
			helix_ids.insert(trees[i][j].second->helix_2_id());
		}

		//if we've not seen this set of helix ids before then add it to our duplicate-less list
		if(helix_id_sets.find(helix_ids) == helix_id_sets.end()){
			helix_id_sets.insert(helix_ids);
			trees_no_dups.push_back(trees[i]);
		}
	}

	TR << "Total trees found: " << trees.size() << std::endl;
	TR << "Total trees without duplicates: " << trees_no_dups.size() << std::endl;

	//***MAKE SURE TREES IS A SUPER-SET OF PATHS***//
	//***CHECK FOR DUPLICATES IN TREE LIST***//

	for(core::Size i=1; i<=trees_no_dups.size(); ++i){

		//if we're finding cycles, repeat the graph
		if(find_cycles){
			ProteinAssemblyGraph::EdgeList cycle=trees_no_dups[i];
			for(core::Size repeat_num=1; repeat_num<repeat_prints; ++repeat_num)
			{
				cycle.insert(cycle.end(), trees_no_dups[i].begin(), trees_no_dups[i].end());
			}

			//If this a 4 helix bundle graph then pop off the last two nodes to get a better repeat
			if(node_size==4)
			{
				cycle.pop_back();
				cycle.pop_back();
			}

			graph->print_tree_simple(cycle);
			if(!simple_print){
				graph->print_tree(cycle, i);
			}
		}

		else{
			TR << "Printing tree " << i << std::endl;
			graph->print_tree_simple(trees_no_dups[i]);
			if(!simple_print){
				graph->print_tree(trees_no_dups[i], i);
			}
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
