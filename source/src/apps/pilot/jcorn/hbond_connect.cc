// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/jcorn
/// @brief  looks at connectedness of h-bonds at an interface
/// @author Jacob Corn (jecorn@u.washington.edu)


#include <devel/init.hh>
#include <core/pose/Pose.hh>

#include <utility/graph/Graph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <protocols/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <basic/MetricValue.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>

#include <protocols/scoring/Interface.hh>
//Utility headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/format.hh>

//C++ headers
#include <iostream>
#include <set>
#include <iterator> // ostream_iterator

//Auto Headers
#include <core/import_pose/import_pose.hh>

//#include <copy>

static basic::Tracer TR( "pilot_apps.interface_graph" );
using namespace core;

void
show_residue_hbonds(
	pose::Pose const & pose,
	Size const seqpos
)
{
	using namespace scoring::hbonds;
	using namespace ObjexxFCL::format;
	HBondSet hbond_set;
	fill_hbond_set( pose, false, hbond_set );

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		HBond const & hb( hbond_set.hbond(i) );
		if ( hb.don_res() == seqpos || hb.acc_res() == seqpos ) {

			std::cout << "RSD_HBOND " <<
				I(4,hb.don_res()) << ' ' << pose.residue( hb.don_res() ).atom_name( hb.don_hatm()) <<
				I(4,hb.acc_res()) << ' ' << pose.residue( hb.acc_res() ).atom_name( hb.acc_atm ()) <<
				F(9,3,hb.energy()) << F(9,3,hb.weight()) << ' ' << hbond_set.allow_hbond(i) << std::endl;
		}
	}

}

void count_residue_hbonds ( scoring::hbonds::HBondSet const & hbond_set, Size const seqpos,  utility::vector1<core::Real> & residue_hbonds ) {
	using namespace scoring::hbonds;

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		HBond const & hb( hbond_set.hbond(i) );
		if ( hb.don_res() == seqpos || hb.acc_res() == seqpos ) {
			residue_hbonds[hb.don_res()] += hb.energy();
			residue_hbonds[hb.acc_res()] += hb.energy();
		}
	}
}

utility::vector1< Real > traverse_for_hbonds ( utility::vector1< utility::vector1< Real > > const & interface_hbond_set ) {
	using namespace core;
	Size const nres( interface_hbond_set.size() );
	TR << "SIZE " << interface_hbond_set.size() << std::endl;
	utility::vector1< Real > link_vector( nres, 0 );

	Real const primary_weight( 2.0 );
	Real const secondary_weight( 1.0 );
	Real const tertiary_weight( 0.5 );

	for ( Size i=1; i<=nres; ++i ) {
		for ( Size j=1; j<=nres; ++j ) {
			if ( i==j ) continue;
			if ( interface_hbond_set[i][j] != 0 ) {
				link_vector[i] += interface_hbond_set[i][j] * primary_weight;

				for ( Size k=1; k<=nres; ++k ) {
					if ( j==k ) continue;
					if ( interface_hbond_set[j][k] != 0 ) {
						link_vector[i] += interface_hbond_set[i][j] * secondary_weight;

						for ( Size l=1; l<=nres; ++l ) {
							if ( k==l ) continue;
							if ( interface_hbond_set[k][l] != 0 ) {
								link_vector[i] += interface_hbond_set[i][j] * tertiary_weight;
							}
						}
					}
				}
			}
		}
	}
	return link_vector;
}

/// @brief writes out a connectivity description of the graph in the dot format
/// readable by things like graphviz. (where the first column "DOT:" should be sed'ed out)
///
/// @param os - [in] - the output stream to write to
void output_interface_graphviz(core::pose::Pose const & pose, utility::graph::Graph const & g, std::ostream & os)
{
	core::Size const jump_num( 1 );
	core::Size const dist_cutoff( 10 );
	protocols::scoring::Interface interface( jump_num );
	interface.distance( dist_cutoff );
	interface.calculate( pose );

	os << "DOT: graph {\n";
	os << "DOT:\tnode[style=filled,color=lightgray];\n";
	os << "DOT:\tgraph[size=8,8];\n";
	// formatting for nodes
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( !interface.is_interface(i) ) continue;
		if ( ! g.get_node(i) ) continue;
		char const name1 = pose.residue( i ).name1();
		core::Size const number = pose.pdb_info()->number( i );
		char const chain_letter = pose.pdb_info()->chain( i );
		core::Size const chain( pose.residue(i).chain() );

		os << "DOT:\t" << name1 << number << "_" << chain_letter << " [ style=filled, color=";
		if ( chain==1 ) os << "lightblue2 ];\n";
		else if ( chain==2 ) os << "goldenrod2 ];\n";
		else os << "darkseagreen ];\n";
	}

	// edge connectivity
	for ( utility::graph::Graph::EdgeListConstIter iter = g.const_edge_list_begin(); iter != g.const_edge_list_end(); ++iter ) {
		core::Size const first_node( (*iter)->get_first_node_ind() );
		core::Size const second_node( (*iter)->get_second_node_ind() );
		if ( first_node == second_node ) continue;
		os << "DOT:\t" << pose.residue( first_node ).name1() << pose.pdb_info()->number( first_node ) << "_" << pose.pdb_info()->chain( first_node );
		os << " -- " << pose.residue( second_node ).name1() << pose.pdb_info()->number( second_node ) << "_" << pose.pdb_info()->chain( second_node ) << ";\n";
	}
	os << "DOT:}\n";
	return;
}


int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		Size jump_num(1);
		Size dist_cutoff( 8 );
		utility::vector1< std::string > pdbnames( basic::options::start_files() );
		for ( utility::vector1< std::string >::const_iterator it = pdbnames.begin(), end = pdbnames.end();
				it != end; ++it ) {

			pose::Pose pose;
			core::import_pose::pose_from_file( pose, *it , core::import_pose::PDB_file);
			core::scoring::ScoreFunctionOP scorefxn =
				core::scoring::get_score_function();
			(*scorefxn)(pose);
			pose.update_residue_neighbors(); // make sure graph_state == GOOD

			protocols::scoring::Interface interface( jump_num );
			interface.distance( dist_cutoff );
			interface.calculate( pose );

			utility::vector1< utility::vector1< Real > > interface_hbond_set;
			for ( Size i=1; i<=pose.size(); ++i ) {
				utility::vector1<Size> zero_vector( pose.size(), 0 );
				interface_hbond_set.push_back( zero_vector );
			}

			/* for debugging
			scoring::hbonds::HBondSet hbond_set1;
			scoring::hbonds::fill_hbond_set( pose, false, hbond_set1 );
			hbond_set1.show(pose);
			*/

			// Pre-minimize sc's to make sure Rosetta is happy and picks up all the hbonds it can
			core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
			mm->clear();
			mm->set_chi( true );
			for ( core::Size i = 1; i <= pose.size(); ++i ) {
				if ( !pose.residue(i).is_protein() ) {
					mm->set_chi( i, false );
					continue;
				}
				// Check for disulfide bonded cysteines
				if ( pose.residue(i).type().is_disulfide_bonded() ) mm->set_chi( i, false );
			}
			protocols::minimization_packing::MinMover min_mover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 1e-5, true/*nblist*/, false/*deriv_check*/  );
			min_mover.apply( pose );

			pose.update_residue_neighbors(); // make sure graph_state == GOOD
			interface.calculate( pose );

			scoring::hbonds::HBondSet hbond_set;
			scoring::hbonds::fill_hbond_set( pose, false, hbond_set );
			// hbond_set2.show(pose);


			utility::graph::Graph g( pose.size() );
			for ( Size i=1; i<=pose.size(); ++i ) {
				if ( interface.is_interface( i ) ) {
					//utility::vector1< Size > i_hbonds( pose.size(), 0 );

					/* std::string calcname("iface_nbrcalc"+utility::to_string(i));
					std::cout << calcname << ": ";
					core::pose::metrics::CalculatorFactory::Instance().register_calculator( calcname, new protocols::pose_metric_calculators::NeighborsByDistanceCalculator( i, dist_cutoff ) );
					basic::MetricValue< std::set<core::Size> > neighbors_mv;
					pose.metric( calcname, "neighbors", neighbors_mv );
					std::ostream_iterator< Size > output( std::cout, " " );
					std::copy( neighbors_mv.value().begin(), neighbors_mv.value().end(), output );
					std::cout << std::endl;
					*/
					count_residue_hbonds( hbond_set, i, interface_hbond_set[i] );
					//TR << i << " " << interface_hbond_set[i] << std::endl;
					for ( core::Size j=1; j<=pose.size(); ++j ) {
						if ( i == j ) continue;
						if ( interface_hbond_set[i][j] > 0 ) {
							if ( !g.get_edge_exists(i, j) ) g.add_edge( i, j );
						}
						if ( interface.is_pair( pose.residue(i), pose.residue(j) ) ) {
							if ( !g.get_edge_exists(i, j) ) g.add_edge( i, j );
						}
					}
				}
			}

			output_interface_graphviz( pose, g, std::cout );


			utility::vector1< Real > link_vector;
			link_vector = traverse_for_hbonds( interface_hbond_set );

			//std::cout << "chain name resi resn nhb" << std::endl;
			for ( Size i=1; i<=pose.size(); ++i ) {
				//if( link_vector[i] < 0 ) {
				std::cout << "HBOND_B: " <<  pose.pdb_info()->name() << " " << pose.pdb_info()->chain( i ) << " " << pose.pdb_info()->number( i ) << " " << pose.residue( i ).name3() << " " << link_vector[ i ] << std::endl;
				//}
			}

		} // for pdbnames
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
