// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DatabaseEntryWorkUnit.cc
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#include <devel/sewing/BundlePairRmsdWorkUnit.hh>


#include <core/pose/util.hh>
#include <protocols/wum/WorkUnitList.hh>
#include <protocols/wum/WorkUnitManager.hh>
#include <protocols/wum/MPI_WorkUnitManager_Slave.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

#include <cstdio>

#include <core/init/init.hh>

//Numeric
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>

// ObjexxFCL libraries
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//TEST
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>

static thread_local basic::Tracer TR( "BundlePairRmsdCalculator" );

namespace BundlePairRmsdCalculator {
	basic::options::IntegerOptionKey const total_fractions( "total_fractions" ); // fraction to run
	basic::options::IntegerOptionKey const fraction( "fraction" ); // fraction to run
	basic::options::IntegerOptionKey const subfraction_size( "subfraction_size" ); // subfraction size
}

struct node_data
{
	int node_index;
	int bundle_index;
	int helix_1_ind;
	int helix_2_ind;
	int n_atoms_in_rms_calc;
};

int
main( int argc, char * argv [] )
{

	try {

	using cppdb::statement;
	using cppdb::result;
	using namespace basic::database::schema_generator;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::wum;

	using namespace core;

	using ObjexxFCL::FArray2D;
	using ObjexxFCL::FArray1D;

	option.add( BundlePairRmsdCalculator::total_fractions, "The number of fractions to split the selection into");
	option.add( BundlePairRmsdCalculator::fraction, "The fraction of results to run RMSD comparisons on");
//	option.add( BundlePairRmsdCalculator::subfraction_size, "Size of each subfraction");

	// initialize core
	core::init::init(argc, argv);

	// Initialize DB
	utility::sql_database::sessionOP db_session(
		basic::database::get_db_session());

	/**Get fraction of comparisons to be run**/
	int fraction = option[BundlePairRmsdCalculator::fraction].def(1);
	int total_fractions = option[BundlePairRmsdCalculator::total_fractions].def(1000);
//	int subfraction_size = option[BundlePairRmsdCalculator::subfraction_size].def(1000000);

	std::string get_count_string =
	"SELECT max(node_id) FROM protein_graph_nodes;";
	cppdb::statement get_count=basic::database::safely_prepare_statement(get_count_string, db_session);
	cppdb::result count_res=basic::database::safely_read_from_database(get_count);

	core::Size total_rows=0;
	while(count_res.next())
	{
		count_res >> total_rows;
	}

	core::Size num_rows_per_fraction = total_rows/total_fractions+1;

	core::Size overall_min_id = num_rows_per_fraction*fraction - num_rows_per_fraction + 1;
	core::Size overall_max_id = num_rows_per_fraction*fraction;
	TR << "overall_min_id: " << overall_min_id << std::endl;
	TR << "overall_max_id: " << overall_max_id << std::endl;

//	core::Size total_sub_fractions=num_rows_per_fraction/subfraction_size+1;
//	core::Size min_id=overall_min_id;
//	core::Size max_id=min_id+subfraction_size;
//
//	TR << "num_rows_per_fraction: " << num_rows_per_fraction << std::endl;
//	TR << "total_sub_fractions: " << total_sub_fractions << std::endl;

	/**Get coordinates**/
	std::string select_helix_coords =
//	"SELECT bh.bundle_id, bh.helix_id, bh.flipped, res.x, res.y, res.z\n"
//	"FROM bundle_helices bh\n"
//	"JOIN structures s ON\n"
//	"	s.struct_id = bh.struct_id\n"
//	"JOIN residue_atom_coords res ON\n"
//	"	bh.struct_id = res.struct_id AND\n"
//	"	res.atomno IN (1,2,3,4) AND\n"
//	"	(res.seqpos BETWEEN bh.residue_begin AND bh.residue_end)\n"
//	"ORDER BY bh.helix_id, res.seqpos, res.atomno;";

	/**For smotifs**/
	"SELECT sp.smotif_pair_id, sss.segment_id, res.x, res.y, res.z\n"
	"FROM smotif_pairs sp\n"
	"JOIN smotifs s1 ON\n"
	"	sp.smotif_id_1 = s1.smotif_id\n"
	"JOIN smotifs s2 ON\n"
	"	sp.smotif_id_2 = s2.smotif_id\n"
	"JOIN secondary_structure_segments sss ON\n"
	"	s1.secondary_struct_segment_id_1 = sss.segment_id OR\n"
	"	s1.secondary_struct_segment_id_2 = sss.segment_id OR\n"
	"	s2.secondary_struct_segment_id_1 = sss.segment_id OR\n"
	"	s2.secondary_struct_segment_id_2 = sss.segment_id\n"
	"JOIN residue_atom_coords res ON\n"
	"	sp.struct_id = res.struct_id AND\n"
	"	res.atomno IN (1,2,3,4) AND\n"
	"	(res.seqpos BETWEEN sss.residue_begin AND sss.residue_end)\n"
	"ORDER BY sss.segment_id, res.seqpos, res.atomno;";

	cppdb::statement coords_stmt=basic::database::safely_prepare_statement(select_helix_coords, db_session);
	cppdb::result coords_res=basic::database::safely_read_from_database(coords_stmt);
	TR << "Done querying coordinates from DB" << std::endl;

	std::map<core::Size, utility::vector1< numeric::xyzVector<core::Real> > > helix_coords;
	std::map<core::Size, utility::vector1< core::Size > > bundle_helices;
	std::map<core::Size, bool > helix_flipped;
	while(coords_res.next())
	{
		core::Size bundle_id, helix_id;
		core::Real x, y, z;

		coords_res >> bundle_id >> helix_id >>  x >> y >> z;
		helix_coords[helix_id].push_back(numeric::xyzVector< core::Real >(x,y,z));
		bundle_helices[bundle_id].push_back(helix_id);
//		helix_flipped[helix_id]=flipped;
	}

	TR << "Done populating helix coords and bundle helices" << std::endl;

	//COUNT NODES
	std::string count_nodes =
		"SELECT\n"
		" count(*)\n"
		"FROM protein_graph_nodes;\n";
	cppdb::statement count_nodes_stmt =
		basic::database::safely_prepare_statement(count_nodes, db_session);
	cppdb::result count_nodes_res =
		basic::database::safely_read_from_database( count_nodes_stmt );

	core::Size n_nodes;
	if ( ! count_nodes_res.next() ) {
		utility_exit_with_message( "Failed to retrieve the number of nodes from the database" );
	}
	count_nodes_res >> n_nodes;
	TR << "Done counting nodes from the database.  n_nodes = " << n_nodes  << std::endl;

	//FILL VECTOR OF NODE DATA
	std::string select_all_nodes =
		"SELECT\n"
		"	node_id, substruct_id, comparison_element_1, comparison_element_2\n"
		"FROM protein_graph_nodes;\n";
	cppdb::statement select_all_nodes_stmt=
		basic::database::safely_prepare_statement(select_all_nodes, db_session);

	cppdb::result select_all_nodes_res=basic::database::safely_read_from_database(select_all_nodes_stmt);

	utility::vector1< node_data > nodes( n_nodes );
	for ( core::Size ii = 1; ii <= n_nodes; ++ii ) {
		if ( ! select_all_nodes_res.next() ) {
			std::cerr << "Retrieved fewer nodes from the select_all_nodes query than the n_nodes query:" << ii << " vs " << n_nodes << std::endl;
			utility_exit_with_message( "problem with BundlePairRMSDCalculatorGPU application" );
		}
		select_all_nodes_res >> nodes[ii].node_index >> nodes[ii].bundle_index >> nodes[ii].helix_1_ind >> nodes[ii].helix_2_ind;
	}

	//INSERT STATEMENT
	std::string comparison_insert =
	"INSERT INTO node_comparisons(node_id_1, node_id_2, rmsd, clash_score) VALUES(?,?,?,?);";
	cppdb::statement comparison_insert_stmt(basic::database::safely_prepare_statement(comparison_insert,db_session));

	//LOOP THROUGH ALL NODES
	db_session->begin();
	for(core::Size ii=1; ii<=n_nodes; ++ii)
	{
		for(core::Size jj=ii+1; jj<=n_nodes; ++jj)
		{
			//ensure the secondary structure elements we are comparing are the right size
			if(helix_coords[ nodes[ii].helix_1_ind ].size() != helix_coords[ nodes[jj].helix_1_ind ].size()) continue;
			if(helix_coords[ nodes[ii].helix_2_ind ].size() != helix_coords[ nodes[jj].helix_2_ind ].size()) continue;

			//ensure we aren't comparing the same smotif to itself. This could be easier if we had the smotif id and not the smotif pair id
			if(helix_coords[ nodes[ii].helix_1_ind ] == helix_coords[ nodes[jj].helix_1_ind ]) continue;
			if(helix_coords[ nodes[ii].helix_2_ind ] == helix_coords[ nodes[jj].helix_2_ind ]) continue;

			utility::vector1< numeric::xyzVector<core::Real> > node_1_coords;
			node_1_coords.insert( node_1_coords.end(),
				helix_coords[ nodes[ii].helix_1_ind ].begin(),
				helix_coords[ nodes[ii].helix_1_ind ].end()
			);
			node_1_coords.insert( node_1_coords.end(),
				helix_coords[ nodes[ii].helix_2_ind ].begin(),
				helix_coords[ nodes[ii].helix_2_ind ].end()
			);

			utility::vector1< numeric::xyzVector<core::Real> > node_2_coords;
			node_2_coords.insert( node_2_coords.end(),
				helix_coords[ nodes[jj].helix_1_ind ].begin(),
				helix_coords[ nodes[jj].helix_1_ind ].end()
			);
			node_2_coords.insert( node_2_coords.end(),
				helix_coords[ nodes[jj].helix_2_ind ].begin(),
				helix_coords[ nodes[jj].helix_2_ind ].end()
			);

			//TESTING
			TR << "Nodes " << nodes[ii].node_index << " " << nodes[jj].node_index << std::endl;
			TR << "Substructure ids " << nodes[ii].bundle_index << " " << nodes[jj].bundle_index << std::endl;
			TR << "Comparison element 1s: " << nodes[ii].helix_1_ind << " " << nodes[jj].helix_1_ind << std::endl;
			TR << "Comparison element 2s: " << nodes[ii].helix_2_ind << " " << nodes[jj].helix_2_ind << std::endl;
//			TR << "helix 1 coords size : " << node_1_coords.size() << std::endl;
//			TR << "helix 2 coords size : " << node_2_coords.size() << std::endl;
			core::Real test_rmsd = numeric::model_quality::calc_rms(node_1_coords,node_2_coords);
			TR << "TEST RMSD: " << test_rmsd << std::endl;

			//Convert to FArrays
			runtime_assert(node_1_coords.size() == node_2_coords.size());

			//Save coords to FArrays
			FArray2D< numeric::Real > p1_coords( 3, node_1_coords.size() );
			FArray2D< numeric::Real > p2_coords( 3, node_2_coords.size() );
			for ( Size i = 1; i <= node_1_coords.size(); ++i )
			{
				for ( Size k = 1; k <= 3; ++k )// k = X, Y and Z
				{
					p1_coords(k,i) = node_1_coords[i](k);
					p2_coords(k,i) = node_2_coords[i](k);
				}
			}

			//uu is rotational matrix that transforms p2coords onto p1coords in a way that minimizes rms. Do this for ONLY the first two helices
			FArray1D< numeric::Real > ww( node_1_coords.size(), 1.0 );//weight matrix, all 1 for my purposes
			FArray2D< numeric::Real > uu( 3, 3, 0.0 );//transformation matrix
			numeric::Real ctx;
			numeric::model_quality::findUU( p1_coords, p2_coords, ww, node_1_coords.size(), uu, ctx );

			//Fill arrays with all coords not involved in the RMSD calculation. These coordinates
			//are used to calculate a clash score.
			utility::vector1< numeric::xyzVector<core::Real> > node_1_other_helix_coords;
//			utility::vector1< core::Size> node_1_helices = bundle_helices[ nodes[ii].bundle_index ];
//			for(core::Size i=1; i<=node_1_helices.size(); ++i)
//			{
//				if( node_1_helices[i] != node_1_helix_id_1 &&
//					node_1_helices[i] != node_1_helix_id_2)
//				{
//					if(helix_flipped[node_1_helices[i]])
//					{
//						node_1_other_helix_coords.insert(node_1_other_helix_coords.end(),
//							helix_coords[node_1_helices[i]].begin(), helix_coords[node_1_helices[i]].end() );
//					}
//					else
//					{
//						utility::vector1< numeric::xyzVector<core::Real> > temp = helix_coords[node_1_helices[i]];
//						reverse(temp.begin(), temp.end());
//						node_1_other_helix_coords.insert( node_1_other_helix_coords.end(), temp.begin(), temp.end() );
//					}
//				}
//			}
//
			utility::vector1< numeric::xyzVector<core::Real> > node_2_other_helix_coords;
//			utility::vector1< core::Size> node_2_helices = bundle_helices[ nodes[jj].bundle_index ];
//			for(core::Size i=1; i<=node_2_helices.size(); ++i)
//			{
//				if( node_2_helices[i] != node_2_helix_id_1 &&
//					node_2_helices[i] != node_2_helix_id_2)
//				{
//					if(helix_flipped[node_2_helices[i]])
//					{
//						node_2_other_helix_coords.insert(node_2_other_helix_coords.end(),
//							helix_coords[node_2_helices[i]].begin(), helix_coords[node_2_helices[i]].end() );
//					}
//					else
//					{
//						utility::vector1< numeric::xyzVector<core::Real> > temp = helix_coords[node_2_helices[i]];
//						reverse(temp.begin(), temp.end());
//						node_2_other_helix_coords.insert( node_2_other_helix_coords.end(), temp.begin(), temp.end() );
//					}
//				}
//			}
//
//			assert(node_1_other_helix_coords.size() == node_2_other_helix_coords.size());

			//Add other helix coords to node vectors, these are now essentiall bundle-vectors
//			node_1_coords.insert( node_1_coords.end(), node_1_other_helix_coords.begin(), node_1_other_helix_coords.end() );
//			node_2_coords.insert( node_2_coords.end(), node_2_other_helix_coords.begin(), node_2_other_helix_coords.end() );

			//Save coords, now with third helix, to new FArrays
			FArray2D< numeric::Real > b1_coords( 3, node_1_coords.size() );
			FArray2D< numeric::Real > b2_coords( 3, node_2_coords.size() );
			for ( Size i = 1; i <= node_1_coords.size(); ++i )
			{
				for ( Size k = 1; k <= 3; ++k )// k = X, Y and Z
				{
					b1_coords(k,i) = node_1_coords[i](k);
					b2_coords(k,i) = node_2_coords[i](k);
				}
			}

			//move bundle to the origin using the center of mass for only the first two helices
			for ( int k = 1; k <= 3; ++k )
			{
				numeric::Real bundle_1_offset = 0.0;
				numeric::Real bundle_2_offset = 0.0;

				for ( int j = 1; j <= node_2_coords.size()-node_1_other_helix_coords.size(); ++j )
				{
					bundle_1_offset += b1_coords(k,j);
					bundle_2_offset += b2_coords(k,j);
				}
				bundle_1_offset /= (node_2_coords.size()-node_1_other_helix_coords.size());
				bundle_2_offset /= (node_2_coords.size()-node_1_other_helix_coords.size());

				for ( int j = 1; j <= node_1_coords.size(); ++j )
				{
					b1_coords(k,j) -= bundle_1_offset;
					b2_coords(k,j) -= bundle_2_offset;
				}
			}

			//transform the coords of the entire bundle using rotational matrix produced by findUU for the first two helices
			FArray2D< numeric::Real > b2_coords_transformed(3, node_2_coords.size());
			for ( int i = 1; i <= node_2_coords.size(); ++i )
			{
				b2_coords_transformed(1,i) = ( uu(1,1)*b2_coords(1,i) )+( uu(1,2)*b2_coords(2,i) ) +( uu(1,3)*b2_coords(3,i) );
				b2_coords_transformed(2,i) = ( uu(2,1)*b2_coords(1,i) )+( uu(2,2)*b2_coords(2,i) ) +( uu(2,3)*b2_coords(3,i) );
				b2_coords_transformed(3,i) = ( uu(3,1)*b2_coords(1,i) )+( uu(3,2)*b2_coords(2,i) ) +( uu(3,3)*b2_coords(3,i) );
			}

			//Calculate RMSD for the first two helices
			numeric::Real tot = 0;
			for ( int i = 1; i <= node_2_coords.size()-node_2_other_helix_coords.size(); ++i )
			{
				for ( int j = 1; j <= 3; ++j )
				{
					tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
				}
			}
			core::Real bundle_pair_rmsd = std::sqrt(tot/(node_2_coords.size()-node_2_other_helix_coords.size()));

//			//calc rms between transformed bundle_2_helix_3 pts and bundle_1_helix_3 pts
//			tot = 0;
//			for(int i = helix_coords[node_1_helix_id_1].size()+helix_coords[node_1_helix_id_2].size()+1;
//				i <= node_2_coords.size(); ++i )
//			{
//				for ( int j = 1; j <= 3; ++j ) {
//					tot += std::pow( b1_coords(j,i) - b2_coords_transformed(j,i), 2 );
//				}
//			}
//			core::Real third_helix_rmsd = std::sqrt(tot/node_1_other_helix_coords.size() );
			core::Real third_helix_rmsd = 0;

			TR << "ALL THE RMSDS: " << bundle_pair_rmsd << " " << test_rmsd << " " << third_helix_rmsd << std::endl;

			comparison_insert_stmt.bind(1, nodes[ii].node_index);
			comparison_insert_stmt.bind(2, nodes[jj].node_index);
			comparison_insert_stmt.bind(3, bundle_pair_rmsd);
			comparison_insert_stmt.bind(4, third_helix_rmsd);
			basic::database::safely_write_to_database(comparison_insert_stmt);
		}
	}
	db_session->commit();
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
