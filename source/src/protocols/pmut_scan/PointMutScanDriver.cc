// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pmut_scan/PointMutScanDriver.cc
/// @brief A protocol that tries to find stability enhancing mutations
/// @author Ron Jacak (ron.jacak@gmail.com)

// Unit headers
#include <protocols/pmut_scan/PointMutScanDriver.hh>
#include <protocols/pmut_scan/Mutant.hh>

//project Headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/graph/Graph.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>

#include <protocols/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/task_operations/RestrictToNeighborhoodOperation.hh>

// Utility Headers
#include <utility>
#include <utility/file/FileName.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

#ifdef USEMPI
/// MPI
#include <mpi.h>
#endif

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace core;
using namespace core::pack::task::operation;
using namespace core::pack::task;

using namespace protocols;
using namespace ObjexxFCL;
using namespace utility;


namespace protocols {
namespace pmut_scan {


static basic::Tracer TR( "protocols.pmut_scan.PointMutScanDriver" );

///
/// @brief
/// Main constructor for the class. What all does it do?
///
PointMutScanDriver::PointMutScanDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string const & list_file, bool output_mutant_structures ) :
	double_mutant_scan_( double_mutant_scan ),
	mutants_list_file_( list_file ),
	output_mutant_structures_( output_mutant_structures ),
	pdb_file_names_( pdb_file_names ),
	DDG_cutoff_(0),
	scorefxn_(core::scoring::get_score_function())
{

#ifdef USEMPI
	tag_ = 1; // need to initialize the tag on all nodes to 1 or MPI_Send/_Recv calls start acting funny
#endif

	int mpi_rank( 0 ), mpi_nprocs( 1 );
#ifdef USEMPI
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );/* get current process id */
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_nprocs );/* get number of processes */
#endif

	MPI_rank_ = (Size)( mpi_rank );
	MPI_nprocs_ = (Size)( mpi_nprocs );

	read_in_structures(); // all processes read in the structures


	// create a scorefxn that will be used for all the mutants
	// (to enable hpatch scoring, the command line weights file flag will have to be used)
	// decompose bb hbond energies into pair energies
	//
	//scoring::ScoreFunctionOP scorefxn = scoring::get_score_function(); //in initialization list
	scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn_->set_energy_method_options( energymethodoptions );

}

///
/// @brief
/// Destructor. What all needs to be done here?
///
PointMutScanDriver::~PointMutScanDriver() {
	//This used to be in the parent application.  For consistency with most JD2-style MPI-compatible apps (which this is not), Finalize has been moved here.
#ifdef USEMPI
	MPI_Finalize();
#endif
}

///
/// @brief
/// Entry point for the pmut_scan protocol.  This is function the app calls to do the scan.
///
void PointMutScanDriver::go() {

	clock_t entire_starttime(0);
	if ( MPI_rank_ == 0 ) {
		// time the protocol, doesn't include the time spent reading in input structures.
		entire_starttime = clock();
	}

	TR << "go(): " << node_name( MPI_rank_ ) << std::endl;

	if ( MPI_rank_ == 0 ) {
		// set up the list of mutations that will be tried.
		// if the user specified a list, then do just those. if not, try all possible combinations of mutants.
		fill_mutations_list();
	}

	barrier(); // do we really want all processes to hold here?
	divide_up_mutations();

	barrier(); // do we really want all processes to hold here?
	make_mutants();

	barrier();
	if ( MPI_rank_ == 0 ) {
		clock_t entire_stoptime = clock();
		TR << "main(): whole protocol took " << ((double)entire_stoptime-entire_starttime) / CLOCKS_PER_SEC << " seconds" << std::endl;
		TR << "go(): DONE with pmut scan." << std::endl;
	}

}

///
std::string PointMutScanDriver::node_name( Size rank ) {

	if ( rank == 0 ) {
		return "master node";
	} else {
		std::stringstream r;
		r << "slave node " << rank;
		return r.str();
	}
}

///
/// Make all processes stop and wait here.
///
void PointMutScanDriver::barrier() {

#ifdef USEMPI
	MPI_Barrier( MPI_COMM_WORLD );
#endif
	std::cout.flush();

}

///
/// @brief
/// Reads in the structure (or list of structures) specified on the command line to a class member variable. Create
/// Pose objects for all of the structures because we'll pass these out to the slave nodes later.
///
/// NOTE: This protocol assumes that if you pass multiple structures, they are all variants of the same structure and you
/// want to use all of them for the ddG calculation.
///
///
void PointMutScanDriver::read_in_structures() {


	// read in all the PDB files into a vector of Pose objects
	//
	utility::vector1< std::string >::iterator input_pdb_filename, last_pdb;
	for ( input_pdb_filename = pdb_file_names_.begin(), last_pdb = pdb_file_names_.end(); input_pdb_filename != last_pdb; ++input_pdb_filename ) {
		pose::Pose pose;
		core::import_pose::pose_from_file( pose, *input_pdb_filename , core::import_pose::PDB_file);
		input_poses_.push_back( pose );
	}

}

///
/// @brief
/// Determines whether the user specified a list of mutants or just wants to do a scan over all possible combinations.
///
/// If we're doing a scan over all possible mutations:
///
/// If we have a two residue protein, 1 and 2, there are 19 possible aa's we can mutate each residue to. Since 1 and 2
/// are independent, the number of possible double mutants is 19*19 = 361. It's easy to enumerate all of the possible
/// double mutants, but we want to come up with an efficient way of making those mutants and calculating the ddGs. The
/// easiest solution is to enumerate all of the possible double mutants, make that a work unit, and then distribute all
/// of the work units out to a large cluster. The downside to this approach is that several processors will end up making
/// the same mutation, at least in part. For example, the double mutant A1C A2C is similar to the mutant A1C A2D. In fact,
/// A1C will have to be paired not only with A2C and A2D, but also A2E, A2F and so on. It would be more efficient to make
/// a pose for A1C, and then go into a second for loop that tries A2C-A2W on that already mutated pose.
///
/// What's the outermost thing this protocol has to do.  For the single mutant scan, you have to try all 19 non-wt aas at
/// every position.  That lends itself to parallelization rather easily.  Each protein position is independent of the others
/// so you can have nres processes each testing the mutations at that position.  At most, each processor will test 19 mutations.
/// With double_mutants, you have to fix one mutation (eg. A1C) and then try all possible other mutations at all other
/// positions.  So, if you have nres processors, each processor will fix some position to 1 of 19 aas and then scan through
/// a mutant at all other positions. Let's assume we have a 10 residue protein.  Position 1 will mutate to 1 of 19 aas.
/// For each of 1 of those 19, we have to test 19 * 9 = 171 other mutations (at the other positions). That results in a
/// grand total of 3249 possibilites.  And that's only residue 1's mutants!  We also have to try to fix the 19 non-wt aas
/// for position 2 and try 19 * 8 = 152 mutations at the other locations for a total of 2888 mutations for just position
/// 2. 3: 19 * 19 * 7 = 2527. 4: 19 * 19 * 6 = 2166.  5: 19 * 19 * 5 = 1805. Continuing on in this fashion leads to a
/// grand grand total of 16245 possible double mutants in a 10 residue protein. Doing the same kind of protocol for a
/// 233 residue protein results in 9,841,221 possible double mutants!
///
/// Testing ~10 million mutants even on 512 cpus could take quite a bit of time. We really need to find a way to prune
/// down the number of possible mutants to just the ones that will be most interesting. I definitely could change it so
/// that if the two mutations are more than some number of Angstroms apart, then don't bother making that mutant and
/// scoring. The question then becomes how often you have a stabilizing first mutant, and then find a stabilizing (better
/// than -0.1) second mutant on the first structure that is more than xAng away. Probably happens often.
///
/// Another problem is that the parallelization is not balanced. Because we have directionality in the approach for
/// testing double mutants - for example, if we've already done 1AC 2AC we don't have to do 2AC 1AC - processor 1 which
/// handles all of the possible mutants at 1 and every other residue has to do way way less
///
///
/// For triple mutants, assuming a 10 residue protein there would be 19 * 19 * 19 * nres(nres+1)/2, or ~377,000, possible
/// mutants. The 233 residue antibody: 186,983,199 possible combinations.
///
///
///
void PointMutScanDriver::fill_mutations_list() {

	if ( !mutants_list_file_.empty() ) {
		read_mutants_list_file( mutants_list_file_ );
		return;
	}

	// otherwise, we're just going to do a scan over all mutations
	// this outer for loop is over all sets of mutations: either single mutants, double mutants, triple mutants, combinations
	// of single, double and triple mutants, etc.
	//utility::vector1< Mutant > stabilizing_mutants;
	//scan_for_mutations( input_poses, scorefxn, stabilizing_mutants, double_mutant_scan_ );

	Size no_double_mutants_possible = 0;
	Size no_double_mutants_excluded_for_distance = 0;
	Size no_double_mutants_excluded_otherwise = 0;
	Size no_single_mutants_excluded_otherwise = 0;


	// use the first structure to determine neighborship for all residues. this neighbor_graph will be used inside the
	// nested for loops to skip mutants that are on opposite sides of the protein.
	utility::vector1< utility::vector1< bool > > neighbors;
	calculate_neighbor_table( input_poses_[1], neighbors );

	pose::Pose & pose = input_poses_[1];
	Size n_residue = pose.size();

	for ( Size resid1 = 1; resid1 <= n_residue; ++resid1 ) {

		// try every type at each position (well, except the native type at this position!)
		for ( Size aa_enum_index_a = 1; aa_enum_index_a <= chemical::num_canonical_aas; ++aa_enum_index_a ) {

			//if ( resid1 > 1 ) { break; } // for debugging only

			if ( pose.residue( resid1 ).aa() == chemical::AA( aa_enum_index_a ) ) { continue; }
			if ( !pose.residue_type( resid1 ).is_protein() ) { continue; }

			MutationData md1(
				pose.residue( resid1 ).name1(),
				oneletter_code_from_aa( chemical::AA( aa_enum_index_a ) ),
				resid1,
				pose.pdb_info()->number( resid1 ),
				pose.pdb_info()->icode( resid1 ),
				pose.pdb_info()->chain( resid1 )
			);

			//single mutant scan
			Mutant m;
			m.add_mutation( md1 ); // the variable mutations is a vector of vectors!
			if ( reject_mutant(m, pose) ) { //offers a chance for child classes to inject mutant selection logic
				++no_single_mutants_excluded_otherwise;
			} else {
				all_mutants_.push_back( m );
				//TR << "fill_mutations_list(): adding mutation: " << m << std::endl;
			}

			// only do a double mutant scan if the user asked for it
			if ( double_mutant_scan_ ) {

				// only need to iterate over higher indexed residues. can't make two mutations at the same position!
				for ( Size resid2 = resid1 + 1; resid2 <= n_residue; ++resid2 ) {

					// check to see if these residues are neighbors of each other. we don't want to make double mutants
					// where the mutants are on opposite sides of the protein.
					if ( neighbors[ resid1 ][ resid2 ] == false ) {
						no_double_mutants_possible += 19;
						no_double_mutants_excluded_for_distance += 19;
						//TR << "skipping residue pair " << md1.mutation_string_PDB_numbering() << " and " << pose.pdb_info()->chain( resid2 ) << "-" << pose.pdb_info()->number( resid2 ) << pose.pdb_info()->icode( resid2 ) << " based on distance" << std::endl;
						continue;
					}

					// try every type at each position (well, except the native type at this position!)
					for ( Size aa_enum_index_b = 1; aa_enum_index_b <= chemical::num_canonical_aas; ++aa_enum_index_b ) {

						if ( pose.residue( resid2 ).aa() == chemical::AA( aa_enum_index_b ) ) { continue; }
						if ( !pose.residue_type( resid2 ).is_protein() ) { continue; }

						no_double_mutants_possible++;

						MutationData md2(
							pose.residue( resid2 ).name1(),
							oneletter_code_from_aa( chemical::AA( aa_enum_index_b ) ),
							resid2,
							pose.pdb_info()->number( resid2 ),
							pose.pdb_info()->icode( resid2 ),
							pose.pdb_info()->chain( resid2 )
						);

						Mutant m;
						m.add_mutation( md1 ); // the variable mutations is a vector of vectors!
						m.add_mutation( md2 ); // the variable mutations is a vector of vectors!
						if ( reject_mutant(m, pose) ) { //offers a chance for child classes to inject mutant selection logic
							++no_double_mutants_excluded_otherwise;
							continue;
						}
						all_mutants_.push_back( m );
						//TR << "fill_mutations_list(): adding mutation: " << m << std::endl;
					}//for all residue types for resid 2
				} // all residues resid2
			}//if a double mutant scan
		}//for all res types for resid 1
	}//for all residues resid1

	if ( MPI_rank_ == 0 ) {
		Size const single_possible = 19 * n_residue;
		TR << "fill_mutations_list(): number single mutants possible: " << single_possible << std::endl;
		TR << "fill_mutations_list(): number single mutants excluded otherwise: " << no_single_mutants_excluded_otherwise << std::endl;
		if ( double_mutant_scan_ ) {
			TR << "fill_mutations_list(): number double mutants possible: " << no_double_mutants_possible << std::endl;
			TR << "fill_mutations_list(): number double mutants excluded for distance: " << no_double_mutants_excluded_for_distance << std::endl;
			TR << "fill_mutations_list(): number double mutants excluded otherwise: " << no_double_mutants_excluded_otherwise << std::endl;
		}
	}

}

///
/// @brief
/// If the user specified mutants, it reads the lines in the mutant list file and parses those lines to get mutation
/// data and then saves them all to the class member variable.
/// Needs access to a pose to translate the lines in the mutations_list file to pose numbering
///
void PointMutScanDriver::read_mutants_list_file( std::string & list_file ) {

	std::ifstream data( list_file.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open mutations file: " + list_file + '\n' );
	}

	// read in all lines in file
	utility::vector1< std::string > mutant_file_lines;
	std::string line;
	while ( getline( data, line ) ) {
		if ( line.size() < 1 || line[0] == '#' ) continue; // skip comment lines
		mutant_file_lines.push_back( line );
	}
	data.close();


	// iterate over all the lines
	for ( Size ii=1; ii <= mutant_file_lines.size(); ++ii ) {
		std::string const & line( mutant_file_lines[ ii ] );
		std::istringstream iss( line );

		char wt_residue, mut_residue, chain;
		std::string position_code;

		Mutant m;

		// there might be more than one mutation per line!
		while ( iss.peek() && !iss.eof() ) {

			iss >> chain >> wt_residue >> position_code >> mut_residue;

			// check to see if an insertion code is present in the position_code string
			// if the string is made of all digits, no icode is present
			Size pdb_resnum; char icode = ' ';
			std::stringstream ss;

			if ( position_code.find_first_not_of("0123456789") == std::string::npos ) {
				icode = ' ';
				ss << position_code;
				ss >> pdb_resnum;

			} else {
				for ( std::string::iterator it = position_code.begin(); it < position_code.end(); ++it ) {
					if ( isdigit(*it) ) {
						ss << (*it);
					} else {
						icode = *it; // assumes that insertion code is only 1-letter!!
					}
				}
				ss >> pdb_resnum; // converts the ss buffer contents to a Size type
			}

			// figure out what the pose residue number for this residue is
			pose::Pose & pose = input_poses_[ 1 ];
			Size pose_resnum = (pose.pdb_info())->pdb2pose( chain, pdb_resnum, icode );

			if ( pose.residue( pose_resnum ).name1() != wt_residue ) {
				TR << "wt_residue: " << wt_residue << ", pdb resnum: " << pdb_resnum << ", pose resnum: " << pose_resnum
					<< ", residue at pose resnum: " << pose.residue( pose_resnum ).name1() << std::endl;
				utility_exit_with_message("Error. Wild-type residue given in mutatons_list file does not match input structure. Please try again.");
			}

			//TR << "Found mutation of " << wt_residue << " to " << mut_residue  << " at position " << pose_resnum << " (pdb chain: '" << chain << "', resnum: '" << pdb_resnum << "', icode: '" << icode << "')" << std::endl;

			MutationData md( wt_residue, mut_residue, pose_resnum, pdb_resnum, icode, chain );
			m.add_mutation( md ); // the variable mutations is a vector of vectors!

		} // done parsing line

		all_mutants_.push_back( m );

	} // end iterating over lines read from input file

}

///
/// @brief
/// Calculates the 10A neighbor graph using the given pose object and then sets values in a 2D array to indicate which
/// resids are neighbors.
///
void PointMutScanDriver::calculate_neighbor_table( pose::Pose & pose, utility::vector1< utility::vector1< bool > > & neighbors ) {

	// size the neighbors 2D table
	neighbors.resize( pose.size(), utility::vector1< bool >( pose.size(), false ) );

	// PointGraph is a one-way graph, which makes it somewhat annoying for iterating over neighbors of a certain
	// position. Only edges to higher-indexed nodes exist. So instead, make a graph which has all the edges at every
	// node to simplify iterating over all neighboring edges.
	core::conformation::PointGraphOP pg( new core::conformation::PointGraph ); // create graph
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg ); // create vertices
	core::conformation::find_neighbors( pg, 10.0 /* Angstrom cutoff */ ); // create edges

	// actually create the neighbor graph from the point graph
	utility::graph::Graph neighbor_graph( pose.size() );
	for ( Size r=1; r <= pose.size(); ++r ) {
		for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex(r).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex(r).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			neighbor_graph.add_edge(r, edge_iter->upper_vertex());
		}
	}

	for ( Size ii=1; ii <= pose.size(); ++ii ) {

		conformation::Residue const & ii_rsd( pose.residue( ii ) );
		for ( utility::graph::EdgeListConstIterator eli = neighbor_graph.get_node( ii )->const_edge_list_begin(),
				eli_end = neighbor_graph.get_node( ii )->const_edge_list_end(); eli != eli_end; ++eli ) {

			Size nb_resnum = (*eli)->get_other_ind( ii );
			if ( nb_resnum < ii ) { continue; } // only want higher indexed residues

			// check to see if any of the atoms on this neighboring residue "interact" with any atoms on the ii residue.
			// our definition of interact: one sc-sc atom pair within 4.5A (BK's suggestion)
			conformation::Residue const & jj_rsd( pose.residue( nb_resnum ) );

			for ( Size jja = jj_rsd.first_sidechain_atom(); jja <= jj_rsd.nheavyatoms(); ++jja ) {
				conformation::Atom const & jja_atom( jj_rsd.atom( jja ) );
				Vector const & jja_atom_xyz = jja_atom.xyz();

				for ( Size iia = ii_rsd.first_sidechain_atom(); iia <= ii_rsd.nheavyatoms(); ++iia ) {
					conformation::Atom const & iia_atom( ii_rsd.atom( iia ) );
					Vector const & iia_atom_xyz = iia_atom.xyz();

					if ( iia_atom_xyz.distance( jja_atom_xyz ) < 4.5 ) {
						neighbors[ ii ][ nb_resnum ] = true; // only set the upper half of the 2D table; i.e. res1 must always be < res2
						break;
					}

				} // ii rsd atoms

				if ( neighbors[ ii ][ nb_resnum ] ) {
					// already found an atom pair within 4.5A; no point in going through all the rest of jj rsd's atoms!
					break;
				}

			} // jj rsd atoms
		}
	}

}

///
/// @brief
/// This function takes the vector of all possible mutants and splits them up as evenly as possible among all the CPUs.
///
void PointMutScanDriver::divide_up_mutations() {

	//TR << "Node " << MPI_rank_ << ", entered method divide_up_mutations()" << std::endl;

	if ( MPI_rank_ == 0 ) {
		//utility::vector1< Mutant > all_mutants_;

		Size const num_mutants_per_cpu = all_mutants_.size() / MPI_nprocs_;
		Size const nextra = all_mutants_.size() - ( num_mutants_per_cpu * MPI_nprocs_ );

		Size my_njobs = ( nextra >= 1 ? 1 : 0 ) + num_mutants_per_cpu;
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			mutants_list_.push_back( all_mutants_[ ii ] );
		}

#ifdef USEMPI
		//TR << "divide_up_mutations(): number of nodes " << MPI_nprocs_ << std::endl;
		Size mutant_offset = my_njobs;

		// send the other nodes their mutations lists so they know what they'll be working on
		for ( Size node_index = 1; node_index < MPI_nprocs_; ++node_index ) {
			Size node_njobs = ( nextra > node_index ? 1 : 0 ) + num_mutants_per_cpu;
			MPI_Send( & node_njobs, 1, MPI_UNSIGNED_LONG, node_index, tag_, MPI_COMM_WORLD );

			for ( Size mutant_index = mutant_offset + 1; mutant_index <= mutant_offset + node_njobs; ++mutant_index ) {
				send_mutant_data_to_node( node_index, all_mutants_[ mutant_index ] );
			}
			mutant_offset += node_njobs;
		}

	} else {
		// slave node. need to receive work order from master node.
		Size my_njobs;
		MPI_Recv( & my_njobs, 1, MPI_UNSIGNED_LONG, 0, tag_, MPI_COMM_WORLD, & stat_ );

		//TR << "divide_up_mutations(): received my_njobs: '" << my_njobs << "'" << std::endl;
		mutants_list_.reserve( my_njobs );
		for ( Size ii = 1; ii <= my_njobs; ++ii ) {
			mutants_list_.push_back( receive_mutant_data_from_node( 0 ) );
		}
#endif
	}

#ifdef USEMPI
	sleep( MPI_rank_ ); // a crude way to order processes...
	for ( Size ii = 1; ii <= mutants_list_.size(); ++ii ) {
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		TR << "divide_up_pdbs(): mutation '" << mutants_list_[ ii ] << "' assigned to " << hostname << " (rank = " << MPI_rank_ << ")" << std::endl;
	}
#endif

}


#ifdef USEMPI
///
///
/// @brief
/// Takes a Mutant and a destination and constructs the MPI_Send call.
///
void PointMutScanDriver::send_mutant_data_to_node( int destination, const protocols::pmut_scan::Mutant & m ) {

	int tag( 1 );

	// each particular mutant can have one, two or more mutations associated with it, make sure to send all of them!
	Size mutant_num_mutations = m.n_mutations();
	//TR << "sending mutant_num_mutations: " << mutant_num_mutations << " to node " << destination << std::endl;
	MPI_Send( & mutant_num_mutations, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );

	for ( utility::vector1< MutationData >::const_iterator iter = m.mutations_begin(); iter != m.mutations_end(); ++iter ) {

		char wt_residue = iter->wt_residue_;
		char mut_residue = iter->mut_residue_;
		//TR << "sending wt_residue: '" << wt_residue << "' and mut_residue: '" << mut_residue << "' to node " << destination << "." << std::endl;
		MPI_Send( & wt_residue, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & mut_residue, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );

		Size pose_resnum = iter->pose_resnum_;
		Size pdb_resnum = iter->pdb_resnum_;
		MPI_Send( & pose_resnum, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & pdb_resnum, 1, MPI_UNSIGNED_LONG, destination, tag, MPI_COMM_WORLD );

		char icode = iter->icode_;
		char chain = iter->chain_;
		//TR << "sending icode: '" << icode << "' and chain: '" << chain << "' to node " << destination << "." << std::endl;
		MPI_Send( & icode, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );
		MPI_Send( & chain, 1, MPI_CHAR, destination, tag, MPI_COMM_WORLD );

	}

}

///
/// @brief
/// Receive mutant data from the master node.  First find out how many mutations are in this mutant and then actually
/// get the mutation data.
///
Mutant PointMutScanDriver::receive_mutant_data_from_node( int source ) {

	int tag( 1 );
	MPI_Status stat;

	Size num_mutations;
	MPI_Recv( & num_mutations, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );
	//TR << "received mutant_num_mutations from node " << source << ": " << num_mutations <<  std::endl;

	Mutant m;
	for ( Size ii = 1; ii <= num_mutations; ++ii ) {

		char wt_residue, mut_residue;
		MPI_Recv( & wt_residue, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & mut_residue, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		//TR << "received wt_residue: " << wt_residue << " and mut_residue: " << mut_residue << " from node " << source << "." << std::endl;

		Size pose_resnum = 1, pdb_resnum = 1;
		MPI_Recv( & pose_resnum, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & pdb_resnum, 1, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD, & stat );

		char icode, chain;
		MPI_Recv( & icode, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		MPI_Recv( & chain, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, & stat );
		//TR << "received icode: '" << icode << "' and chain: '" << chain << "' from node " << source << "." << std::endl;

		MutationData md( wt_residue, mut_residue, pose_resnum, pdb_resnum, icode, chain );
		m.add_mutation( md );
	}

	//TR << "receive_mutant_data_from_node(): received mutant '" << m << "'" << std::endl;

	return m;
}

#endif

///
/// @brief
/// Calls make_specific_mutant on all mutants assigned to this node.
/// Also responsible for creating the score function that's used for all mutants.
///
void PointMutScanDriver::make_mutants() {

	utility::vector1< pose::Pose > mutant_poses( input_poses_.size() ); // this will get set in the function below
	utility::vector1< pose::Pose > native_poses( input_poses_.size() );

	// print out a header to the terminal
	if ( MPI_rank_ == 0 ) {
		TR << format::A( "mutation" ) << format::X(3) << format::A( "mutation_PDB_numbering" ) << format::X(3) << format::A( "average_ddG" ) << format::X(3) << format::A( "average_total_energy" ) << std::endl;
	}

	for ( Size ii=1; ii <= mutants_list_.size(); ++ii ) {
		Mutant & m = mutants_list_[ ii ];

		//TR << "make_mutants(): making mutant: " << m << std::endl;

		// the make specific mutant function changes both the mutant and native poses. we want to start with our
		// original starting structures each time though. so we have to copy the input poses to some working native
		// and mutant poses vectors.
		for ( Size ii=1; ii <= input_poses_.size(); ++ii ) {
			mutant_poses[ ii ] = input_poses_[ ii ];
			native_poses[ ii ] = input_poses_[ ii ];
		}

		make_specific_mutant( mutant_poses, native_poses, m, "", "" );
		// this will result in the Mutant object 'm' being modified, and since m is a reference, the original mutants_list_
		// will be modified, as well.
	}

}

///
/// @brief
/// Function which takes in a mutation to make (either single, double or more) and calls itself recursively until the desired
/// structure is created. Useful for testing certain combinations of mutations (like putting two double mutants together)
/// without having to run an entire scan protocol that would cover those mutations.
///
void PointMutScanDriver::make_specific_mutant( utility::vector1< pose::Pose > & mutant_poses, utility::vector1< pose::Pose > & native_poses,
	Mutant & m, std::string mutation_string, std::string mutation_string_PDB_numbering ) {

	//TR << "make_specific_mutant() called. mutant_poses.size(): " << mutant_poses.size() << ", native_poses.size(): " << native_poses.size()
	// << ", num mutations: " << m.n_mutations() << ", mutation_string: " << mutation_string << std::endl;

	// if the mutants vector has more than element, we have to take out the first element of the vector
	if ( m.n_mutations() > 1 ) {

		// need to make the first mutation and call this function recursively
		MutationData md = m.pop_mutation();

		// make the first mutation on the mutant_poses
		for ( Size ii = 1; ii <= native_poses.size(); ++ii ) {

			// make the specific mutation, but don't do any scoring; the scorefxn is needed for packing
			make_mutant_structure( mutant_poses[ ii ], native_poses[ ii ], md );

		}

		std::stringstream out;
		out << mutation_string;
		if ( mutation_string != "" ) { out << ","; }
		out << md.mutation_string();
		std::string updated_mutation_string = out.str();
		out.str("");

		out << mutation_string_PDB_numbering;
		if ( mutation_string_PDB_numbering != "" ) { out << ","; }
		out << md.mutation_string_PDB_numbering();
		std::string updated_mutation_string_PDB_numbering = out.str();


		make_specific_mutant( mutant_poses, native_poses, m, updated_mutation_string, updated_mutation_string_PDB_numbering );

	} else {
		// make the last mutation, calculate the ddG, and print out the results
		MutationData md = m.pop_mutation();

		//TR << "make_specific_mutant(): making final mutation: " << md << std::endl;

		Energy sum_mutant_scores = 0.0;
		Energy average_mutant_score = 0.0;

		Energy sum_native_scores = 0.0;
		Energy average_native_score = 0.0;

		utility::vector1< Real > native_poses_total_energies( native_poses.size() );
		utility::vector1< Real > mutant_poses_total_energies( native_poses.size() );

		for ( Size ii=1; ii <= native_poses.size(); ++ii ) {
			// make the specific mutation, but don't do any scoring; the scorefxn is needed for packing
			// send in the input_pose for the mutant. that way the mutant poses will be "returned" because mutant_poses
			// is actually a reference!
			make_mutant_structure( mutant_poses[ii], native_poses[ii], md );

			// score the created mutant structure
			pose::Pose & mutant_pose = mutant_poses[ ii ];
			Energy mutant_score = score( mutant_pose );
			mutant_poses_total_energies[ ii ] = mutant_score;
			sum_mutant_scores += mutant_score;

			// score the update native structure
			pose::Pose & native_pose = native_poses[ ii ];
			Energy native_score = score( native_pose );
			native_poses_total_energies[ ii ] = native_score;
			sum_native_scores += native_score;

			if ( output_mutant_structures_ ) {
				std::stringstream out;
				out << mutation_string;
				if ( mutation_string != "" ) { out << ", "; }
				out << md.mutation_string();
				utility::file::FileName fn( pdb_file_names_[ ii ] );
				std::string mutant_filename = fn.base() + "." + out.str() + ".pdb";
				mutant_pose.dump_scored_pdb( mutant_filename, *scorefxn_ );
			}

		}

		average_mutant_score = sum_mutant_scores / mutant_poses.size();
		average_native_score = sum_native_scores / native_poses.size();

		Real ddG_mutation = average_mutant_score - average_native_score;
		if ( ddG_mutation > DDG_cutoff_ ) {
			return;
		}

		std::stringstream out;
		out << mutation_string;
		if ( mutation_string != "" ) { out << ","; }
		out << md.mutation_string();
		std::string final_mutation_string = out.str();

		out.str("");
		out << mutation_string_PDB_numbering;
		if ( mutation_string_PDB_numbering != "" ) { out << ","; }
		out << md.mutation_string_PDB_numbering();
		std::string final_mutation_string_PDB_numbering = out.str();


		TR << final_mutation_string << format::X(3) << final_mutation_string_PDB_numbering << format::X(3) << format::F( 9,3,ddG_mutation ) << format::X(3) << format::F( 9,2,average_mutant_score ) << std::endl;


		/*TR << "native poses total energies: ";
		for ( Size ii=1; ii <= native_poses_total_energies.size(); ++ii ) {
		TR << native_poses_total_energies[ ii ] << ", ";
		}
		TR << std::endl;
		TR << "mutant poses total energies: ";
		for ( Size ii=1; ii <= mutant_poses_total_energies.size(); ++ii ) {
		TR << mutant_poses_total_energies[ ii ] << ", ";
		}
		TR << std::endl;*/
		TR.flush_all_channels();


	} // end loop over all mutants

}

///
/// @brief
/// Given mutant and native pose references and the mutation to make, this function constructs all the necessary PackerTask
/// Operations and Movers to apply the mutation and repacking steps to both the mutant and native poses.
///
void PointMutScanDriver::make_mutant_structure( pose::Pose & mutant_pose, pose::Pose & native_pose, MutationData const & md ) {

	Size resid = md.pose_resnum();
	chemical::AA mut_aa = chemical::aa_from_oneletter_code( md.mut_residue() );

	// need to create a neighborhood by distance calculator so we can identify neighbors of the mutated residue
	std::stringstream out;
	out << md.mutation_string() << "_mutant_nb_calculator";
	std::string calculator_name = out.str();

	pose::metrics::PoseMetricCalculatorOP mutant_nb_calculator( new pose_metric_calculators::NeighborsByDistanceCalculator( resid ) );
	pose::metrics::CalculatorFactory::Instance().register_calculator( calculator_name, mutant_nb_calculator );

	basic::MetricValue< std::set< Size > > mv_neighbors;
	mutant_pose.metric( calculator_name, "neighbors", mv_neighbors );
	std::set< Size > const neighbor_set( mv_neighbors.value() );

	//TR << "make_mutant_structure(): neighbor_set: ";
	//for ( std::set< Size >::iterator it = neighbor_set.begin() ; it != neighbor_set.end(); it++ ) {
	// TR << *it << ", ";
	//}
	//TR << std::endl;

	TaskFactoryOP native_tf( new TaskFactory() );
	TaskFactoryOP mutant_tf( new TaskFactory() );

	// the restrict operation class (which in the end is just a TaskOperation) takes a calculator during construction. I've already
	// created that calculator above.  This operation will disable repacking and design at all positions except those in the neighborhood
	// of the mutated position.
	TaskOperationCOP nb_op( new task_operations::RestrictToNeighborhoodOperation( calculator_name ) );
	native_tf->push_back( nb_op ); mutant_tf->push_back( nb_op );

	// extra task operations we want to also include
	// the restrict residue to repacking ops are used to make sure that only repacking and not design is done to the residues in the neighborhood
	InitializeFromCommandlineOP init_op( new InitializeFromCommandline() );
	native_tf->push_back( init_op ); mutant_tf->push_back( init_op );

	IncludeCurrentOP ic_op( new IncludeCurrent() );
	native_tf->push_back( ic_op ); mutant_tf->push_back( ic_op );

	RestrictResidueToRepackingOP mutant_repack_op( new RestrictResidueToRepacking() );
	RestrictResidueToRepackingOP wt_repack_op( new RestrictResidueToRepacking() ); // will include one extra residue to repack
	for ( Size ii = 1; ii <= mutant_pose.size(); ++ii ) {
		// resid is the position on the original pose. ii is the position on the copy.
		if ( ii == resid ) {
			// do design on this position
			utility::vector1< bool > keep_canonical_aas( chemical::num_canonical_aas, false );
			keep_canonical_aas[ mut_aa ] = true;
			RestrictAbsentCanonicalAASOP restrict_op( new RestrictAbsentCanonicalAAS( ii, keep_canonical_aas ) );
			mutant_tf->push_back( restrict_op );
			wt_repack_op->include_residue( ii ); // for the wild type, don't design on the mutant resid - but do allow repacking
		} else {
			// make this position repackable only; because of the commutativity of packer task ops, only the residues that are in the neighborhood
			// of the mutant will be allowed to repack. the restrict to neighborhood op will disallow packing at all positions not near the mutant.
			mutant_repack_op->include_residue( ii );
			wt_repack_op->include_residue( ii );
		}
	}
	native_tf->push_back( wt_repack_op );
	mutant_tf->push_back( mutant_repack_op );

	//TR << "Finished creating all TaskOperation's and TaskFactory's. Creating MoveMap." << std::endl;

	kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	std::set< core::Size >::const_iterator iter, end;
	for ( iter = neighbor_set.begin(), end = neighbor_set.end(); iter != end; ++iter ) {
		//movemap_->set_bb(i, true); // don't do any backbone minimization
		movemap->set_chi( *iter, true ); // but do minimize the side chains
	}
	//movemap->show( std::cout, mutant_pose.size() );

	//TR << "Movemap created... Beginning repacking/minimization of mutant pose." << std::endl;

	// create an actual PackerTask from the TaskFactory
	pack::task::PackerTaskOP scan_task = mutant_tf->create_task_and_apply_taskoperations( mutant_pose );
	//scan_task->num_to_be_packed();
	//TR << "mutant packer task: " << *scan_task << std::endl;  // generates a TON of output

	// now create the movers that will do the repacking and minimization
	protocols::minimization_packing::PackRotamersMoverOP mutant_repacker_mover( new protocols::minimization_packing::PackRotamersMover( scorefxn_, scan_task, 2 ) ); // ndruns: 2
	protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover( movemap, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true ) ); // use nb_list: true
	protocols::minimization_packing::TaskAwareMinMoverOP task_aware_min_mover( new protocols::minimization_packing::TaskAwareMinMover( min_mover, mutant_tf ) );
	protocols::moves::SequenceMoverOP seq_mover( new protocols::moves::SequenceMover );
	seq_mover->add_mover( mutant_repacker_mover );
	seq_mover->add_mover( task_aware_min_mover );

	seq_mover->apply( mutant_pose );

	//TR << "Beginning repacking/minimization of wt pose." << std::endl;

	// create an actual PackerTask from the TaskFactory
	pack::task::PackerTaskOP wt_task = native_tf->create_task_and_apply_taskoperations( native_pose );

	// now create the movers that will do the repacking and minimization of the native structure
	protocols::minimization_packing::PackRotamersMoverOP native_pack_mover( new protocols::minimization_packing::PackRotamersMover( scorefxn_, wt_task, 2 ) ); // ndruns: 2
	min_mover = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( movemap, scorefxn_, option[ OptionKeys::run::min_type ].value(), 0.01, true ) ); // use nb_list: true
	task_aware_min_mover = protocols::minimization_packing::TaskAwareMinMoverOP( new protocols::minimization_packing::TaskAwareMinMover( min_mover, native_tf ) );
	seq_mover = protocols::moves::SequenceMoverOP( new protocols::moves::SequenceMover );
	seq_mover->add_mover( native_pack_mover );
	seq_mover->add_mover( task_aware_min_mover );

	seq_mover->apply( native_pose );

	// this needs to get recreated each time around
	pose::metrics::CalculatorFactory::Instance().remove_calculator( calculator_name );
	mutant_nb_calculator = nullptr;

	return;

} // done with make_mutant_structure

//Virtual functions, refactored out so they can be overridden by child AlterSpecDisruptionDriver

/// @brief score the pose for the purposes of determining if a mutation is "good" or not.  In the base implementation, it's just a scorefunction call, but in child implementations it may be fancier (for example, calculating a binding energy instead)
core::Energy PointMutScanDriver::score(core::pose::Pose & pose) {
	return (*scorefxn_)(pose);
}

// setters used by the unit tests only
void PointMutScanDriver::set_ddG_cutoff( Real threshold ) {
	DDG_cutoff_ = threshold;
}

///
/// @brief
/// returns a const iterator to the beginning of the Mutant data member variable vector
///
utility::vector1< Mutant >::const_iterator PointMutScanDriver::mutants_begin() const {
	return all_mutants_.begin();
}

///
/// @brief
/// returns a const iterator to the end of the Mutant data member variable vector
///
utility::vector1< Mutant >::const_iterator PointMutScanDriver::mutants_end() const {
	return all_mutants_.end();
}

///
/// @brief
/// returns the size of the Mutant data member variable vector
///
Size PointMutScanDriver::n_mutants() const {
	return all_mutants_.size();
}

core::scoring::ScoreFunctionCOP PointMutScanDriver::get_scorefxn() const {
	return scorefxn_;
}


} // namespace pmut_scan
} // namespace protocols
