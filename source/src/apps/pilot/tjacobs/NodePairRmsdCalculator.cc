// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

//Basic
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

//Core
#include <devel/init.hh>
#include <core/types.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

// ObjexxFCL libraries
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

//C++
#include <cstdio>

static basic::Tracer TR( "NodePairRmsdCalculator" );

struct node_data
{
	core::Size node_id;
	core::Size substructure_id; //Can be a bundle_id, sandwich_id, or smotif_pair id
	core::Size comparison_element_1; //Can be a helix_id, a strand_id, or an smotif_id
	core::Size comparison_element_2; //Can be a helix_id, a strand_id, or an smotif_id
	// core::Size n_atoms_in_rms_calc;
};

//Map from element_id to coordinates
typedef std::map< core::Size, utility::vector1<numeric::xyzVector<core::Real> > > ElementCoords;


//Create a vector of the data for all nodes
utility::vector1<node_data>
get_nodes(
	utility::sql_database::sessionOP db_session
){
	using cppdb::statement;
	using cppdb::result;

	std::string count_nodes =
		//BUNDLES AND SMOTIFS
		//  "SELECT\n"
		//  " count(*)\n"
		//  "FROM protein_graph_nodes;\n";

		//STRANDS
		"SELECT\n"
		"\tcount(*)\n"
		"FROM strand_pairs_graph;";

	statement count_nodes_stmt =
		basic::database::safely_prepare_statement(count_nodes, db_session);
	result count_nodes_res =
		basic::database::safely_read_from_database( count_nodes_stmt );

	if ( ! count_nodes_res.next() ) {
		utility_exit_with_message( "Failed to retrieve the number of nodes from the database" );
	}
	core::Size n_nodes;
	count_nodes_res >> n_nodes;
	TR << "Done counting nodes from the database.  n_nodes = " << n_nodes  << std::endl;

	std::string select_all_nodes =
		//BUNDLES AND SMOTIFS
		//  "SELECT\n"
		//  " node_id, substruct_id, comparison_element_1, comparison_element_2\n"
		//  "FROM protein_graph_nodes;\n";

		//STRANDS
		"SELECT\n"
		"\tstrand_pairs_graph_id as node_id,\n"
		"\tsandwich_graph_id as substructure_id,\n"
		"\tstrand_graph_near_id_1 as comparison_element_1,\n"
		"\tstrand_graph_near_id_2 as comparison_element_2\n"
		"FROM strand_pairs_graph;";

	statement select_all_nodes_stmt=
		basic::database::safely_prepare_statement(select_all_nodes, db_session);

	result select_all_nodes_res=basic::database::safely_read_from_database(select_all_nodes_stmt);


	utility::vector1< node_data > nodes( n_nodes );
	for ( core::Size ii = 1; ii <= n_nodes; ++ii ) {
		if ( ! select_all_nodes_res.next() ) {
			std::cerr << "Retrieved fewer nodes from the select_all_nodes query than the n_nodes query:" << ii << " vs " << n_nodes << std::endl;
			utility_exit_with_message( "problem with BundlePairRMSDCalculatorGPU application" );
		}
		select_all_nodes_res >>
			nodes[ii].node_id >>
			nodes[ii].substructure_id >>
			nodes[ii].comparison_element_1 >>
			nodes[ii].comparison_element_2;
	}

	return nodes;
}

//Create a mapping from substructure_id to a list of element_ids
std::map< core::Size, utility::vector1<core::Size> >
get_substructure_map(
	utility::sql_database::sessionOP db_session
){
	using cppdb::statement;
	using cppdb::result;

	std::string select_substructures =
		//BUNDLES AND SMOTIFS

		//STRANDS
		"SELECT\n"
		"\tswg.sandwich_graph_id as substructure_id,\n"
		"\tsg.strand_graph_id as element_id\n"
		"FROM sandwich_graph swg\n"
		"JOIN strand_pairs_graph spg ON\n"
		"\tswg.strand_pairs_graph_id_i = spg.strand_pairs_graph_id OR\n"
		"\tswg.strand_pairs_graph_id_j = spg.strand_pairs_graph_id\n"
		"JOIN strand_graph sg ON\n"
		"\tspg.strand_graph_near_id_1 = sg.strand_graph_id OR\n"
		"\tspg.strand_graph_near_id_2 = sg.strand_graph_id;";

	statement select_substructures_stmt=
		basic::database::safely_prepare_statement(select_substructures, db_session);

	result select_all_nodes_res=basic::database::safely_read_from_database(select_substructures_stmt);

	std::map< core::Size, utility::vector1<core::Size> > substructure_map;
	while ( select_all_nodes_res.next() )
			{
		core::Size substructure_id, element_id;
		select_all_nodes_res >> substructure_id >> element_id;
		substructure_map[substructure_id].push_back(element_id);
	}

	return substructure_map;
}

//Create a mapping from element_id to list of coordinates
ElementCoords
get_element_coords(
	utility::sql_database::sessionOP db_session
)
{
	std::string select_element_coords =
		//BUNDLES
		// "SELECT bh.bundle_id, bh.helix_id, res.x, res.y, res.z\n"
		// "FROM bundle_helices bh\n"
		// "JOIN structures s ON\n"
		// " s.struct_id = bh.struct_id\n"
		// "JOIN residue_atom_coords res ON\n"
		// " bh.struct_id = res.struct_id AND\n"
		// " res.atomno IN (1,2,3,4) AND\n"
		// " (res.seqpos BETWEEN bh.residue_begin AND bh.residue_end)\n"
		// "ORDER BY bh.helix_id, res.seqpos, res.atomno;";

		//SMOTIFS
		// "SELECT sp.smotif_pair_id, sss.segment_id, res.x, res.y, res.z\n"
		// "FROM smotif_pairs sp\n"
		// "JOIN smotifs s1 ON\n"
		// " sp.smotif_id_1 = s1.smotif_id\n"
		// "JOIN smotifs s2 ON\n"
		// " sp.smotif_id_2 = s2.smotif_id\n"
		// "JOIN secondary_structure_segments sss ON\n"
		// " s1.secondary_struct_segment_id_1 = sss.segment_id OR\n"
		// " s1.secondary_struct_segment_id_2 = sss.segment_id OR\n"
		// " s2.secondary_struct_segment_id_1 = sss.segment_id OR\n"
		// " s2.secondary_struct_segment_id_2 = sss.segment_id\n"
		// "JOIN residue_atom_coords res ON\n"
		// " sp.struct_id = res.struct_id AND\n"
		// " res.atomno IN (1,2,3,4) AND\n"
		// " (res.seqpos BETWEEN sss.residue_begin AND sss.residue_end)\n"
		// "ORDER BY sss.segment_id, res.seqpos, res.atomno;";

		//STRANDS
		"SELECT sg.strand_graph_id element_id, res.x, res.y, res.z\n"
		"FROM strand_graph sg\n"
		"JOIN structures s ON\n"
		"\ts.struct_id = sg.struct_id\n"
		"JOIN residue_atom_coords res ON\n"
		"\ts.struct_id = res.struct_id AND\n"
		"\tres.atomno IN (1,2,3,4) AND\n"
		"\t(res.seqpos BETWEEN sg.residue_start AND sg.residue_end)\n"
		"ORDER BY sg.strand_graph_id, res.seqpos, res.atomno;";

	cppdb::statement coords_stmt=basic::database::safely_prepare_statement(select_element_coords, db_session);
	cppdb::result coords_res=basic::database::safely_read_from_database(coords_stmt);
	TR << "Done querying coordinates from DB" << std::endl;

	ElementCoords element_coords;
	while ( coords_res.next() )
			{
		core::Size element_id;
		core::Real x, y, z;

		coords_res >> element_id >>  x >> y >> z;
		element_coords[element_id].push_back(numeric::xyzVector< core::Real >(x,y,z));
	}

	return element_coords;
}

std::pair<core::Real, core::Real>
calc_rmsd_and_clash_score(
	node_data const & node_1,
	node_data const & node_2,
	ElementCoords & element_coords,
	std::map<core::Size, utility::vector1<core::Size> > & substructure_map
){
	using namespace ObjexxFCL;
	assert(element_coords[node_1.comparison_element_1].size() == element_coords[node_2.comparison_element_1].size());
	assert(element_coords[node_1.comparison_element_2].size() == element_coords[node_2.comparison_element_2].size());

	core::Size element_1_size = element_coords[node_1.comparison_element_1].size();
	core::Size element_2_size = element_coords[node_1.comparison_element_2].size();
	core::Size num_rmsd_atoms = element_1_size + element_2_size;

	utility::vector1<core::Size> const & substructure_1_elements =
		substructure_map[node_1.substructure_id];

	core::Size substructure_1_natoms=0;
	for ( utility::vector1<core::Size>::const_iterator element_it=substructure_1_elements.begin();
			element_it != substructure_1_elements.end(); ++element_it ) {
		substructure_1_natoms+=element_coords[(*element_it)].size();
	}

	utility::vector1<core::Size> const & substructure_2_elements =
		substructure_map[node_2.substructure_id];

	core::Size substructure_2_natoms=0;
	for ( utility::vector1<core::Size>::const_iterator element_it=substructure_2_elements.begin();
			element_it != substructure_2_elements.end(); ++element_it ) {
		substructure_2_natoms+=element_coords[(*element_it)].size();
	}

	//Save coords to FArrays
	FArray2D< numeric::Real > p1_coords( 3, substructure_1_natoms);
	FArray2D< numeric::Real > p2_coords( 3, substructure_2_natoms);

	core::Size substructure_1_coord_counter=1;
	core::Size substructure_2_coord_counter=1;
	for ( core::Size ii = 1; ii <= element_1_size; ++ii ) {
		for ( core::Size k = 1; k <= 3; ++k ) {// k = X, Y and Z
			p1_coords(k,substructure_1_coord_counter) = element_coords[node_1.comparison_element_1][ii](k);
			p2_coords(k,substructure_2_coord_counter) = element_coords[node_2.comparison_element_1][ii](k);
		}
		++substructure_1_coord_counter;
		++substructure_2_coord_counter;
	}

	for ( core::Size ii = 1; ii <= element_2_size; ++ii ) {
		for ( core::Size k = 1; k <= 3; ++k ) {// k = X, Y and Z
			p1_coords(k,substructure_1_coord_counter) = element_coords[node_1.comparison_element_2][ii](k);
			p2_coords(k,substructure_2_coord_counter) = element_coords[node_2.comparison_element_2][ii](k);
		}
		++substructure_1_coord_counter;
		++substructure_2_coord_counter;
	}

	//Add the rest of the coordinates from substructure 1 to p1_coords
	for ( utility::vector1<core::Size>::const_iterator element_it=substructure_1_elements.begin();
			element_it != substructure_1_elements.end(); ++element_it ) {
		if ( (*element_it) != node_1.comparison_element_1 && (*element_it) != node_1.comparison_element_2 ) {
			for ( core::Size ii = 1; ii <= element_coords[(*element_it)].size() ; ++ii ) {
				for ( core::Size k = 1; k <= 3; ++k ) {// k = X, Y and Z
					p1_coords(k,substructure_1_coord_counter) = element_coords[(*element_it)][ii](k);
				}
				++substructure_1_coord_counter;
			}
		}
	}

	//Add the rest of the coordinates from substructure 2 to p2_coords
	for ( utility::vector1<core::Size>::const_iterator element_it=substructure_2_elements.begin();
			element_it != substructure_2_elements.end(); ++element_it ) {
		if ( (*element_it) != node_2.comparison_element_1 && (*element_it) != node_2.comparison_element_2 ) {
			for ( core::Size ii = 1; ii <= element_coords[(*element_it)].size() ; ++ii ) {
				for ( core::Size k = 1; k <= 3; ++k ) {// k = X, Y and Z
					p2_coords(k,substructure_2_coord_counter) = element_coords[(*element_it)][ii](k);
				}
				++substructure_2_coord_counter;
			}
		}
	}

	//uu is rotational matrix that transforms p2coords onto p1coords in a way that minimizes rms.
	//Do this for the two element-pairs in nodes i and j
	FArray1D< numeric::Real > ww( num_rmsd_atoms, 1.0 );//weight matrix, all 1 for my purposes
	FArray2D< numeric::Real > uu( 3, 3, 0.0 );//Rotational matrix
	numeric::Real ctx;
	numeric::model_quality::findUU( p1_coords, p2_coords, ww, num_rmsd_atoms, uu, ctx );

	//move substructure to the origin using the center of mass for only the first two helices
	for ( core::Size k = 1; k <= 3; ++k ) {
		numeric::Real structure_1_com = 0.0;
		numeric::Real structure_2_com = 0.0;

		for ( core::Size j = 1; j <= num_rmsd_atoms; ++j ) {
			structure_1_com += p1_coords(k,j);
			structure_2_com += p2_coords(k,j);
		}
		structure_1_com /= num_rmsd_atoms;
		structure_2_com /= num_rmsd_atoms;

		for ( core::Size j = 1; j <= substructure_1_natoms; ++j ) {
			p1_coords(k,j) -= structure_1_com;
		}

		for ( core::Size j = 1; j <= substructure_2_natoms; ++j ) {
			p2_coords(k,j) -= structure_2_com;
		}
	}

	//Transform p2_coords using the UU matrix defined above
	// FArray2D< numeric::Real > p2_coords_transformed(3, node_2_coords.size());
	for ( core::Size i = 1; i <= substructure_1_natoms; ++i ) {
		p2_coords(1,i) = ( uu(1,1)*p2_coords(1,i) )+( uu(1,2)*p2_coords(2,i) ) +( uu(1,3)*p2_coords(3,i) );
		p2_coords(2,i) = ( uu(2,1)*p2_coords(1,i) )+( uu(2,2)*p2_coords(2,i) ) +( uu(2,3)*p2_coords(3,i) );
		p2_coords(3,i) = ( uu(3,1)*p2_coords(1,i) )+( uu(3,2)*p2_coords(2,i) ) +( uu(3,3)*p2_coords(3,i) );
	}

	//Calculate RMSD for the two element pairs
	numeric::Real tot = 0;
	for ( core::Size k = 1; k <= 3; ++k ) {
		for ( core::Size i = 1; i <= num_rmsd_atoms; ++i ) {
			tot += std::pow( p1_coords(k,i) - p2_coords(k,i), 2 );
		}
	}
	core::Real element_pair_rmsd = std::sqrt(tot/(num_rmsd_atoms));

	//Calculate clash score for all elements not in RMSD comparison
	core::Real min_dist_sq=1000000;
	core::Real dist_sq;
	for ( core::Size i=num_rmsd_atoms+1; i<=substructure_1_natoms; ++i ) {
		for ( core::Size j=num_rmsd_atoms+1; j<=substructure_2_natoms; ++j ) {
			dist_sq=0;
			for ( int k = 1; k <= 3; ++k ) {
				dist_sq+=std::pow( p1_coords(k,i) - p2_coords(k,j), 2 );
			}
			min_dist_sq = std::min(min_dist_sq, dist_sq);
		}
	}

	return std::make_pair(element_pair_rmsd, min_dist_sq);
}

int
main( int argc, char * argv [] )
{
	using utility::vector1;
	using namespace ObjexxFCL;

	// initialize core
	devel::init(argc, argv);

	// Initialize DB
	utility::sql_database::sessionOP db_session(
		basic::database::get_db_session());

	vector1<node_data> nodes =
		get_nodes(db_session);

	std::map<core::Size, vector1<core::Size> > substructure_map =
		get_substructure_map(db_session);

	ElementCoords element_coords =
		get_element_coords(db_session);

	std::string comparison_insert =
		"INSERT INTO node_comparisons(node_id_1, node_id_2, rmsd, clash_score) VALUES(?,?,?,?);";
	cppdb::statement comparison_insert_stmt(basic::database::safely_prepare_statement(comparison_insert,db_session));
	db_session->begin();
	for ( core::Size i=1; i<=nodes.size(); ++i ) {
		for ( core::Size j=i+1; j<=nodes.size(); ++j ) {
			//Don't compare nodes from the same substructure
			if ( nodes[i].substructure_id != nodes[j].substructure_id ) {
				std::pair<core::Real, core::Real> rmsd_and_clash =
					calc_rmsd_and_clash_score(nodes[i], nodes[j], element_coords, substructure_map);
				comparison_insert_stmt.bind(1, nodes[i].node_id);
				comparison_insert_stmt.bind(2, nodes[j].node_id);
				comparison_insert_stmt.bind(3, rmsd_and_clash.first);
				comparison_insert_stmt.bind(4, rmsd_and_clash.second);
				basic::database::safely_write_to_database(comparison_insert_stmt);
			}
		}
	}
	db_session->commit();

	return 0;
}
