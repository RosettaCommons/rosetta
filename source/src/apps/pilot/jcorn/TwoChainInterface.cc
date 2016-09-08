// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FlexDesign
/// @brief jcorn


/* example usage:

*/

// libRosetta headers


#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

// // C++ headers
#include <cstdlib>


//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;

//using utility::vector1;

static THREAD_LOCAL basic::Tracer TT( "pilot_apps.jcorn.two_chain_interface", basic::t_trace );
static THREAD_LOCAL basic::Tracer TD( "pilot_apps.jcorn.two_chain_interface", basic::t_debug );
static THREAD_LOCAL basic::Tracer TI( "pilot_apps.jcorn.two_chain_interface", basic::t_info );

class TwoChainInterface {

public:
	TwoChainInterface();	//  empty constructor
	TwoChainInterface( pose::Pose &, Size const, Size const ); // constructor with chain numbers
	TwoChainInterface( pose::Pose &, char const, char const ); // constructor with chain letters

	// returns a std::pair of vector1s of residues on one side of interface. pair->first[x] interacts with pair->second[x]
	std::pair< utility::vector1<Size>, utility::vector1<Size> >  interface( )  { return list_interface_; }

private:
	Size chain1_number_, chain2_number_;
	char chain1_letter_, chain2_letter_;
	Size ch1_begin_, ch1_end_;
	Size ch2_begin_, ch2_end_;
	Size partner1_res_, partner2_res_;
	utility::vector1< Size > chain1_interface_, chain2_interface_;
	utility::vector1< Size >::iterator size_iter_;
	std::pair < utility::vector1<Size>, utility::vector1<Size> > list_interface_;

	void get_chain_terminii_( pose::Pose const & pose, Size const chain1_id, Size const chain2_id)	// figures out chain terminii from chain numbers
	{
		chain1_number_ = chain1_id;
		chain2_number_ = chain2_id;
		ch1_begin_ = pose.conformation().chain_begin( chain1_number_ );
		ch1_end_ = pose.conformation().chain_end( chain1_number_ );
		ch2_begin_ = pose.conformation().chain_begin( chain2_number_ );
		ch2_end_ = pose.conformation().chain_end( chain2_number_ );
	}

	Size chain_letter_to_number ( pose::Pose const & pose, char const chain_id )	// helper method to convert a single chain letter to a single chain number. returns 0 if chain letter not found
	{
		char temp_letter_ = chain_id;
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.pdb_chains()[i] == temp_letter_ ) {
				return pose.chain( i );
			}
			else continue;
		}
		return 0;	// if this happens, we haven't found our chain letter in the pose
	}


	void find_interface( pose::Pose & pose, char const chain1_id, char const chain2_id )	// overloading to take chain letters
	{
		chain1_letter_ = chain1_id;
		chain2_letter_ = chain2_id;

		// convert everything in search_chain_letters to rosetta-style chain numbers
		chain1_number_ = chain_letter_to_number( pose, chain1_letter_ );
		chain2_number_ = chain_letter_to_number( pose, chain2_letter_ );

		find_interface( pose, chain1_number_, chain2_number_ );	// Here's where the magic happens
	}


	void find_interface( pose::Pose & pose, Size const chain1_id, Size const chain2_id )		// MAJOR INTERFACE FINDING HAPPENS HERE!
	{
		chain1_number_ = chain1_id;
		chain2_number_ = chain2_id;
		get_chain_terminii_( pose, chain1_number_, chain2_number_ );

		pose.update_residue_neighbors();	// make sure graph_state == GOOD
		scoring::TenANeighborGraph const & tenA_neighbor_graph = pose.energies().tenA_neighbor_graph();

		for ( partner1_res_ = ch1_begin_; partner1_res_ <= ch1_end_; ++partner1_res_ ) {  // for all nodes in chain1 == for all residues in chain 1
			for ( graph::Graph::EdgeListConstIter edge_iter = tenA_neighbor_graph.get_node( partner1_res_ )->const_upper_edge_list_begin();
				  edge_iter != tenA_neighbor_graph.get_node( partner1_res_ )->const_upper_edge_list_end();
				  ++edge_iter ) {   // for all edges of node
				partner2_res_ = (*edge_iter)->get_other_ind( partner1_res_ );	// get node on other edge of that node == 2nd residue index
				if ( ( partner2_res_ >= ch2_begin_ ) && (partner2_res_ <= ch2_end_ ) ) {   // if that node(residue) is in chain 2
					chain1_interface_.push_back( partner1_res_ );	// add partner1 residue
					chain2_interface_.push_back( partner2_res_ );	// add partner2 residue
				}
				else continue;
			} // END for all edges of node
		}	// END for all nodes in chain1
		list_interface_ = std::make_pair( chain1_interface_, chain2_interface_ );	// populate list_interface
	}	// END find_interface( pose, chain1_id, chain2_id ) method

};		// END CLASS

TwoChainInterface::TwoChainInterface( pose::Pose & pose, Size const chain1_id, Size const chain2_id ) {	// constructor with arguments
	chain1_number_ = chain1_id;
	chain2_number_ = chain2_id;
	find_interface( pose, chain1_number_, chain2_number_ );
}

TwoChainInterface::TwoChainInterface( pose::Pose & pose, char const chain1_id, char const chain2_id ) {	// constructor with arguments
	chain1_letter_ = chain1_id;
	chain2_letter_ = chain2_id;
	find_interface( pose, chain1_letter_, chain2_letter_ );
}

TwoChainInterface::TwoChainInterface() { }	// empty constructor

int
main( int argc, char* argv[] )
{
	try {

	pose::Pose pose;
	devel::init( argc, argv );

	utility::vector1 < char > search_chain_letters;
	utility::vector1 < Size > search_chain_numbers;
	TwoChainInterface two_chain_interface;
	std::pair< utility::vector1<Size>, utility::vector1<Size> > interface;

	core::import_pose::pose_from_file( pose, "/Users/jcorn/svn/workspaces/jcorn/LJ111.pdb" , core::import_pose::PDB_file);

	//utility::vector1 < char > pdb_chains = pose.pdb_chains();
	search_chain_letters.push_back('A');
	search_chain_letters.push_back('B');
	search_chain_letters.push_back('C');

	// method to search every chain in search_chain_letters against every other chain in search_chain_letters
	for (utility::vector1< char >::iterator chain1 = search_chain_letters.begin(); chain1 != search_chain_letters.end(); ++chain1) {
		for (utility::vector1< char >::iterator chain2 = chain1+1; chain2 != search_chain_letters.end(); ++chain2) {
			two_chain_interface = TwoChainInterface( pose, *chain1, *chain2 );
			interface = two_chain_interface.interface();
			TI << *chain1 << *chain2 << " " << interface.first << std::endl;
			TI << *chain1 << *chain2 << " " << interface.second << std::endl;
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}	// END main

