// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_multithreaded_ig_correctness.cc
/// @brief A cxx11thread application that tests that multithreaded interaction graph computation yields an
/// identical interaction graph to that produced by the single-threaded interaction graph computation.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>

// protocol headers
#include <protocols/symmetry/SetupForSymmetryMover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR("apps.pilot.vmullig.test_multithreaded_ig_correctness");

OPT_KEY (Integer, expected_node_count)
OPT_KEY (Integer, expected_edge_count)

/// @brief Indicate which options are relevant.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( multithreading::total_threads );
	option.add_relevant( multithreading::interaction_graph_threads );
	option.add_relevant( packing::repack_only );
	option.add_relevant( in::file::s );
	option.add_relevant( symmetry::symmetry_definition );
	NEW_OPT( expected_node_count, "Expected number of nodes.  Not used if not provided.", 0);
	NEW_OPT( expected_edge_count, "Expected number of edges.  Not used if not provided..", 0);
}

/// @brief Are two packer energy values equal to within a (hard-coded) cutoff?
bool
values_equal( core::PackerEnergy const val1, core::PackerEnergy const val2 ) {
	return std::abs( val1 - val2 ) < 1.0e-4;
}

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		TR << "Starting test_multithreaded_ig_correctness." << std::endl;
		TR << "Pilot application created 27 July 2019 by Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org)." << std::endl;
		TR << "This application is intended for unit and integration testing, not for production." << std::endl;

		static const std::string errmsg( "Error in test_multithreaded_ig_correctness application: " );

#ifdef MULTI_THREADED

		// Get options
		if ( ! (option [ in::file::s ].user() && option[ in::file::s ]().size() == 1 ) ) {
			utility_exit_with_message( errmsg + "Please specify exactly one input PDB file with the -in:file:s option.");
		}
		std::string const filename( option[ in::file::s ]()[1] );
		core::Size const nthreads( option[multithreading::total_threads]() );
		core::Size const npackthreads( option[multithreading::interaction_graph_threads]() );
		runtime_assert_string_msg( nthreads > 1, errmsg + "The total number of threads specified with the -multithreading:total_threads option must be greater than 1." );
		runtime_assert_string_msg( npackthreads > 1, errmsg + "The number of packing threads specified with the -multithreading:interaction_graph_threads option must be greater than 1." );
		if ( option[ expected_node_count ].user() ) {
			runtime_assert_string_msg( option[ expected_node_count ]() > 0, errmsg + "The expected node count must be positive." );
		}
		if ( option[ expected_edge_count ].user() ) {
			runtime_assert_string_msg( option[ expected_edge_count ]() > 0, errmsg + "The expected edge count must be positive." );
		}
		core::Size const expected_nodes( option[ expected_node_count ]() ), expected_edges( option[ expected_edge_count ]() );

		// Construct the pose that we'll be packing.
		TR << "Importing pose from " << filename << "." << std::endl;
		core::pose::PoseOP pose( core::import_pose::pose_from_file( filename ) );
		runtime_assert_string_msg( pose->total_residue() > 0, errmsg + "Import of pose from " + filename + " failed." );

		// If a symmetry definition is specified, symmetrize the pose.
		if ( option[ symmetry::symmetry_definition ].user() ) {
			std::string const symmfile( option[ symmetry::symmetry_definition ].value() );
			TR << "Setting up symmetry from file " << symmfile << "." << std::endl;
			protocols::symmetry::SetupForSymmetryMover setsymm( symmfile );
			setsymm.apply( *pose );
		}

		// Delete the following -- for debugging only.
		//pose->dump_pdb("testpose.pdb");

		//Create the scorefunction and the packretask:
		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function( true ) );
		core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( *pose, 1 /*Request one thread.*/ ) );
		core::pack::task::PackerTaskOP task_multithreaded( core::pack::task::TaskFactory::create_packer_task( *pose, npackthreads ) ); //Request npackthreads threads, or nthreads if nthreads is less.
		if ( option[ packing::repack_only ].value() ) {
			task->restrict_to_repacking();
			task_multithreaded->restrict_to_repacking();
		}

		//Generate rotamers:
		core::pack::rotamer_set::RotamerSetsOP rotsets( core::pack::rotamer_set::RotamerSetsFactory::create_rotamer_sets( *pose ) );

		// Compute the interaction graph with one thread.
		TR << "Computing interaction graph with a single thread." << std::endl;
		core::pose::PoseOP pose_copy1( pose->clone() );
		core::pack::interaction_graph::AnnealableGraphBaseOP intgraph1(nullptr);
		core::pack::pack_rotamers_setup( *pose_copy1, *sfxn, task, rotsets, intgraph1 );
		core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pc_intgraph1( utility::pointer::dynamic_pointer_cast<core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph>(intgraph1) );
		runtime_assert_string_msg( pc_intgraph1 != nullptr, errmsg + "The interaction graph computed on a single thread is not a PrecomputedPairEnergiesInteractionGraph." );

		// Compute the interaction graph with N threads.
		TR << "Computing interaction graph with " << npackthreads << " threads." << std::endl;
		core::pose::PoseOP pose_copy2( pose->clone() );
		core::pack::interaction_graph::AnnealableGraphBaseOP intgraph2(nullptr);
		core::pack::pack_rotamers_setup( *pose_copy2, *sfxn, task_multithreaded, rotsets, intgraph2 );
		core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pc_intgraph2( utility::pointer::dynamic_pointer_cast<core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph>(intgraph2) );
		runtime_assert_string_msg( pc_intgraph2 != nullptr, errmsg + "The interaction graph computed on a single thread is not a PrecomputedPairEnergiesInteractionGraph." );

		// Check that the two graphs match.
		core::Size edgecount1(0), edgecount2(0);
		utility::vector1< std::string > problems;
		core::Size const nodecount1( intgraph1->get_num_nodes() ), nodecount2( intgraph2->get_num_nodes() );
		if ( expected_nodes != 0 ) {
			if ( nodecount1 != expected_nodes ) {
				problems.push_back( std::string("The single-threaded interaction graph had " + std::to_string(nodecount1) + " nodes, but we expected " + std::to_string(expected_nodes) + "." )  );
			} else {
				TR << "The single-threaded interaction graph has the expected number of nodes." << std::endl;
			}
			if ( nodecount2 != expected_nodes ) {
				problems.push_back( std::string("The multi-threaded interaction graph had " + std::to_string(nodecount2) + " nodes, but we expected " + std::to_string(expected_nodes) + "." )  );
			} else {
				TR << "The single-threaded interaction graph has the expected number of nodes." << std::endl;
			}
		}
		if ( nodecount1 != nodecount2 ) {
			problems.push_back( std::string( "The single-threaded interaction graph has " + std::to_string( nodecount1 ) + "nodes, while the multi-threaded interaction graph has " + std::to_string( nodecount2 ) + "." ) );
		} else /*( nodecount1 == nodecount2 )*/ {
			TR << "Checking " << nodecount1  << " interaction graph nodes." << std::endl;
			for ( core::Size i(1); i<=nodecount1; ++i ) {
				core::Size const numstates1( intgraph1->get_num_states_for_node(i) ), numstates2( intgraph2->get_num_states_for_node(i) );
				if ( numstates1 != numstates2 ) {
					problems.push_back( std::string( "On node " + std::to_string(i) + ", the single-threaded interaction graph has " + std::to_string( numstates1 ) + " states, but the multi-threaded interaction graph has " + std::to_string( numstates2 ) + "." ) );
				} else {
					// Check onebody energies match:
					core::pack::interaction_graph::PrecomputedPairEnergiesNode const * node1(
						dynamic_cast< core::pack::interaction_graph::PrecomputedPairEnergiesNode const * >(pc_intgraph1->get_node(i))
					);
					core::pack::interaction_graph::PrecomputedPairEnergiesNode const * node2(
						dynamic_cast< core::pack::interaction_graph::PrecomputedPairEnergiesNode const * >(pc_intgraph2->get_node(i))
					);
					if ( ( node1 != nullptr ) && ( node2 != nullptr ) ) {
						for ( core::Size istate(1); istate <= numstates1; ++istate ) {
							core::PackerEnergy const energy1( node1->get_one_body_energy(istate) );
							core::PackerEnergy const energy2( node2->get_one_body_energy(istate) );
							if ( !values_equal( energy1, energy2 ) ) {
								problems.push_back( std::string( "For node " + std::to_string(i) + ", state " + std::to_string( istate ) + ", the single-threaded onebody energy is " + std::to_string( energy1 ) + ", while the multi-threaded onebody energy is " + std::to_string( energy2 ) + "." ) );
							}
						}
					} else {
						if ( node1 == nullptr ) {
							problems.push_back( std::string( "Node " + std::to_string(i) + " of the single-threaded interaction graph could not be cast to a PrecomputedPairEnergiesNode!" ) );
						}
						if ( node2 == nullptr ) {
							problems.push_back( std::string( "Node " + std::to_string(i) + " of the multi-threaded interaction graph could not be cast to a PrecomputedPairEnergiesNode!" ) );
						}
					}
				}

				for ( core::Size j(i+1); j<=nodecount1; ++j ) {
					bool const has_edge1( pc_intgraph1->get_edge_exists( i, j ) );
					bool const has_edge2( pc_intgraph2->get_edge_exists( i, j ) );
					if ( has_edge1 ) ++edgecount1;
					if ( has_edge2 ) ++edgecount2;
					if ( has_edge1 != has_edge2 ) {
						if ( has_edge1 ) {
							problems.push_back( std::string( "The single-threaded interaction graph has an edge between nodes " + std::to_string(i) + " and " + std::to_string(j) + " but the multithreaded graph does not." ) );
						} else {
							problems.push_back( std::string( "The multithreaded interaction graph has an edge between nodes " + std::to_string(i) + " and " + std::to_string(j) + " but the single-threaded graph does not." ) );
						}
					} else {
						if ( has_edge1 && has_edge2 ) {
							core::pack::interaction_graph::PrecomputedPairEnergiesEdge const * edge1(
								dynamic_cast< core::pack::interaction_graph::PrecomputedPairEnergiesEdge const * >( pc_intgraph1->find_edge( i, j ) )
							);
							core::pack::interaction_graph::PrecomputedPairEnergiesEdge const * edge2(
								dynamic_cast< core::pack::interaction_graph::PrecomputedPairEnergiesEdge const * >( pc_intgraph2->find_edge( i, j ) )
							);
							runtime_assert_string_msg( edge1 != nullptr, errmsg + "Single-threaded interaction graph edge between nodes " + std::to_string(i) + " and " + std::to_string(j) + " is not a PrecomputedPairEnergiesEdge!" );
							runtime_assert_string_msg( edge2 != nullptr, errmsg + "Multi-threaded interaction graph edge between nodes " + std::to_string(i) + " and " + std::to_string(j) + " is not a PrecomputedPairEnergiesEdge!" );
							for ( core::Size istate(1), istatemax(intgraph1->get_num_states_for_node(i)); istate<=istatemax; ++istate ) {
								for ( core::Size jstate(1), jstatemax(intgraph1->get_num_states_for_node(j)); jstate<=jstatemax; ++jstate ) {
									core::PackerEnergy const energy1( edge1->get_two_body_energy(istate,jstate) );
									core::PackerEnergy const energy2( edge2->get_two_body_energy(istate,jstate) );
									if ( !values_equal( energy1, energy2 ) ) {
										problems.push_back( std::string( "The single-threaded graph edge between node " + std::to_string(i) + ", state " + std::to_string(istate) + " and node " + std::to_string(j) + ", state " + std::to_string(jstate) + " has energy " + std::to_string(energy1) + ", while the multi-threaded equivalent has energy " + std::to_string(energy2) + "."  ) );
									}
								}
							}
						}
					}
				}
			}
		}

		//Check that the total edge count matched:
		if ( edgecount1 != edgecount2 ) {
			problems.push_back ( std::string( "The single-threaded interaction graph had " + std::to_string(edgecount1) + " edges, but the multi-threaded graph had " + std::to_string( edgecount2 ) + " edges." ) );
		} else {
			TR << "Checked " << edgecount1 << " edges." << std::endl;
		}

		//Check that, if the user has specified expected edge counts, the actual edge counts match:
		if ( expected_edges != 0 ) {
			if ( edgecount1 != expected_edges ) {
				problems.push_back( std::string("The single-threaded interaction graph had " + std::to_string(edgecount1) + " edges, but we expected " + std::to_string(expected_edges) + "." )  );
			} else {
				TR << "The single-threaded interaction graph has the expected number of edges." << std::endl;
			}
			if ( edgecount2 != expected_edges ) {
				problems.push_back( std::string("The multi-threaded interaction graph had " + std::to_string(edgecount2) + " edges, but we expected " + std::to_string(expected_edges) + "." )  );
			} else {
				TR << "The multi-threaded interaction graph has the expected number of edges." << std::endl;
			}
		}

		// Throw an error if there are problems:
		if ( problems.size() > 0 ) {
			std::string errmsg_full( errmsg + "The following errors were found:");
			for ( core::Size i(1), imax(problems.size()); i<=imax; ++i ) {
				errmsg_full += "\n" + problems[i];
			}
			errmsg_full += "\n\nExiting with error status!\n";
			utility_exit_with_message( errmsg_full );
		} else {
			TR << "No problems found!" << std::endl;
		}

#else // not MULTI_THREADED

		utility_exit_with_message( errmsg + "The test_multithreaded_ig_correctness application can only be run in the multithreaded build of Rosetta.  Please compile with the extras=cxx11thread option." );

#endif //MULTI_THREADED

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	TR << "Application test_multithreaded_ig_correctness completed successfully.  Exiting with status 0 (no errors)." << std::endl;
	TR.flush();

	return 0;
}
