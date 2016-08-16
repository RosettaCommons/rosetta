// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// libRosetta headers
//#include <basic/options/option.hh>

#include <core/types.hh>


#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <devel/init.hh>

#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>


#include <protocols/optimize_weights/OptEData.fwd.hh>
#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/OptEMultifunc.hh>

#include <core/graph/Graph.hh>


#include <core/io/pdb/pdb_writer.hh>


#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>


#include <protocols/simple_moves/PackRotamersMover.hh>


#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/file/file_sys_util.hh>
//#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>


//silly using/typedef

using namespace core;
using namespace scoring;
using namespace optimization;
using namespace protocols::optimize_weights;

using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


utility::vector1< std::string >
get_native_pdb_names();


void
simple_opte_test();

void
optimize_weights(
	utility::vector1< std::string > const & filenames,
	utility::vector1< std::string > const & native_filenames,
	core::scoring::EnergyMap const & free_parameters,
	core::scoring::EnergyMap const & fixed_parameters,
	utility::vector1< Real > const & starting_reference_energies,
	core::scoring::EnergyMap & optimized_free_parameters,
	utility::vector1< Real > & optimized_reference_energies
);

void
iterative_optE();


Real
measure_sequence_recovery(
	utility::vector1< std::string > const & native_pdb_names,
	utility::vector1< std::string > const & names_for_output_pdbs,
	core::scoring::ScoreFunctionOP sfxn
);

Size
opte_num_inner_iterations(
	Size outer_loop_counter
);

Real
opte_weight_mixing_factor(
	Size outer_loop_counter,
	Size inner_loop_counter
);

core::scoring::ScoreFunctionOP
create_score_function_from_weights_and_refEs
(
	core::scoring::EnergyMap const & emap,
	utility::vector1< Real > const & reference_energies
);

void
initialize_free_and_fixed(
	core::scoring::EnergyMap & free_parameters,
	core::scoring::EnergyMap & fixed_parameters
);

bool
converged(
	core::scoring::EnergyMap & free_parameters_prev,
	core::scoring::EnergyMap & free_parameters_curr,
	utility::vector1< Real > const & reference_energies_prev,
	utility::vector1< Real > const & reference_energies_curr
);


void
setup_pdbnames_next_round(
	Size const outer_loop_counter,
	utility::vector1< std::string  > & pdbs_next_round,
	utility::vector1< std::string > const & native_pdb_names
);

bool
accept_new_weight_set(
	Size outer_loop_counter,
	Size inner_loop_counter,
	Real latest_seq_recovery_rate,
	Real last_sequence_recovery_rate,
	Real inner_loop_scale_factor
);

void
write_parameters_to_std_out(
	core::scoring::EnergyMap & free_parameters,
	utility::vector1< Real > const & reference_energies
);

utility::vector1< std::string >
get_native_pdb_names()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// read list file of pdbs
	utility::vector1< std::string > filenames;


	std::string const listfile( start_file() );
	std::ifstream data( listfile.c_str() );
	std::string line;
	while ( getline( data,line ) ){
		filenames.push_back( line );
	}
	data.close();

	return filenames;
}


void
get_opte_data(
	pose::Pose & pose,
	pose::Pose & native_pose,
	ScoreFunction const & scorefxn,
	ScoreTypes & score_list,
	ScoreTypes & fixed_score_vec,
	OptEData & opte_data
)
{
	using namespace pose;
	using namespace pack::rotamer_set;

	pack::task::PackerTaskOP dummy_task = pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->set_bump_check( false );
	dummy_task->or_include_current( true );
	dummy_task->temporarily_fix_everything();

	scorefxn.setup_for_packing( pose, dummy_task->repacking_residues(), dummy_task->designing_residues() );

	utility::vector1< bool > task_mask( pose.total_residue(), false );
	Size num_diffs_betwee_native_and_input( 0 );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );
		task_mask[ ii ] = true;
		task->restrict_to_residues( task_mask );
		task_mask[ ii ] = false;

		PNatAAOptEPositionDataOP this_pos_data = new PNatAAOptEPositionData;
		this_pos_data->set_position( ii );
		this_pos_data->set_native_aa( native_pose.residue( ii ).aa() );
		if ( native_pose.residue( ii ).aa() != pose.residue( ii ).aa() ) {
			//std::cout << "native_residue # " << ii << " of " << native_pose.residue( ii ).aa() << " differs with pose residue " << pose.residue( ii ).aa() << std::endl;
			++num_diffs_betwee_native_and_input;
		}
		this_pos_data->set_neighbor_count(
					pose.energies().tenA_neighbor_graph().get_node( ii )->num_neighbors_counting_self() );

		graph::GraphCOP packer_neighbor_graph(  pack::create_packer_graph( pose, scorefxn, task ) );

		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( ii ) );
//		RotamerSetOP rotset = RotamerSetFactory::create_rotamer_set( pose.residue( ii ) );

		rotset->set_resid( ii );
		rotset->build_rotamers( pose, scorefxn, *task, packer_neighbor_graph );
		scorefxn.prepare_rotamers_for_packing( pose, *rotset );

		// First, need a vector of energy maps
		utility::vector1< EnergyMap > emap_vector( rotset->num_rotamers() );

		// Call the new energy map fn
		rotset->compute_one_body_energy_maps( pose, scorefxn, *task, packer_neighbor_graph, emap_vector );

		for ( Size jj = 1; jj <= rotset->num_rotamers(); ++jj ) {

			EnergyMap & emap_total( emap_vector[jj] );

			// Hacky limit for fa_rep
			if( emap_total[ fa_rep ] > 10.0 ) emap_total[ fa_rep ] = 10.0;

			utility::vector1< Real > energy_info;
			utility::vector1< Real > fixed_energy_info;

			for( utility::vector1< ScoreType >::iterator score_type_iter = score_list.begin(),
							end_iter = score_list.end() ; score_type_iter != end_iter ; ++score_type_iter ) {
				energy_info.push_back( emap_total[ *score_type_iter ] );
			}

			for( utility::vector1< ScoreType >::iterator score_type_iter = fixed_score_vec.begin(),
							end_iter = fixed_score_vec.end() ; score_type_iter != end_iter ; ++score_type_iter ) {
				fixed_energy_info.push_back( emap_total[ *score_type_iter ] );
			}

			PNatAAOptERotamerDataOP new_rot_line = new PNatAAOptERotamerData(
				(*rotset->rotamer( jj )).aa(), jj,
				energy_info, fixed_energy_info );

			this_pos_data->add_rotamer_line_data( new_rot_line );
		}
		// Done with rotamers for this position, store this position data object
		opte_data.add_position_data( this_pos_data );
	}
	std::cout << "num_diffs_betwee_native_and_input: " << num_diffs_betwee_native_and_input << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
simple_opte_test()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace io::pdb;


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace optimization;
	using namespace pack::rotamer_set;

	ScoreFunction scorefxn;

	/*
	scorefxn.set_weight( fa_atr, 1.00 );
	scorefxn.set_weight( fa_rep, 1.00 );
	scorefxn.set_weight( fa_sol, 1.00 );
	scorefxn.set_weight( fa_dun, 1.00 );
	scorefxn.set_weight( fa_pair, 1.00 );
	scorefxn.set_weight( p_aa_pp, 1.00 );
	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );
	*/

	EnergyMap rand_emap;
	rand_emap[ vdw ] = numeric::random::rg().uniform();
	rand_emap[ env ] = numeric::random::rg().uniform();
	rand_emap[ pair ] = numeric::random::rg().uniform();
	rand_emap[ cbeta ] = numeric::random::rg().uniform();
	rand_emap[ rama ] = numeric::random::rg().uniform();
	rand_emap[ p_aa_pp ] = numeric::random::rg().uniform();
	rand_emap[ cenpack ] = numeric::random::rg().uniform();

	scorefxn.set_weight( vdw,     rand_emap[ vdw ] );
	scorefxn.set_weight( env,     rand_emap[ env ] );
	scorefxn.set_weight( pair,    0.33  );
	scorefxn.set_weight( cbeta,   rand_emap[ cbeta ]  );
	scorefxn.set_weight( rama,    rand_emap[ rama ]  );
	scorefxn.set_weight( p_aa_pp, rand_emap[ p_aa_pp ] );
	scorefxn.set_weight( cenpack, rand_emap[ cenpack ] );
	scorefxn.set_weight( rg,      3.0 );
	scorefxn.set_weight( hs_pair, 1.0 );
	scorefxn.set_weight( ss_pair, 1.0 );
	scorefxn.set_weight( rsigma,  1.0 );
	scorefxn.set_weight( sheet,   1.0 );


	// Make two EnergyMaps - one is the set of weights
	// which are included in the optimization, the other
	// indicates which should be fixed at a certain non-zero
	// value;

	EnergyMap include_terms( scorefxn.weights() );
	EnergyMap fixed_terms;

	/*
	fixed_terms.set( fa_atr, 0.8 );
	*/

	fixed_terms.set( pair,     0.33 );
	fixed_terms.set( rg,      3.0 );
	//fixed_terms.set( cenpack, 1.0 );
	fixed_terms.set( hs_pair, 1.0 );
	fixed_terms.set( ss_pair, 1.0 );
	fixed_terms.set( rsigma,  1.0 );
	fixed_terms.set( sheet,   1.0 );

	// Loop through the energy maps and figure out how
	// many terms are included, and how many total are free.

	int include_count( 0 );
	int fixed_count( 0 );
	int score_index( 1 );
	for( EnergyMap::const_iterator include_itr = include_terms.begin(),
			end_include_itr = include_terms.end(),
			fixed_itr = fixed_terms.begin() ;
			include_itr != end_include_itr ; ++include_itr, ++fixed_itr, ++score_index ) {
		if( (*include_itr) != 0.0 ) {
			++include_count;
		}
		if( (*fixed_itr) != 0.0 ) {
			++fixed_count;
			if( (*include_itr) == 0.0 ) {
				utility_exit_with_message( "ERROR:  Fixed energy term not included!" );
			}
		}
	}

	int free_count( include_count - fixed_count );

	std::cout << "Including " << include_count << " energy terms" << std::endl;
	std::cout << "Fixing " << fixed_count << " energy terms" << std::endl;
	std::cout << "Optimizing " << free_count << " energy terms" << std::endl;

	ScoreTypes score_list;
	for( int i=1 ; i <= n_score_types ; ++i ) {
		if( include_terms[ ScoreType(i) ] != 0.0 ) {
			score_list.push_back( ScoreType(i) );
		}
	}

	ScoreTypes fixed_score_list;
	for( int i=1 ; i <= n_score_types ; ++i ) {
		if( fixed_terms[ ScoreType(i) ] != 0.0 ) {
			fixed_score_list.push_back( ScoreType(i) );
		}
	}

	for( ScoreTypes::const_iterator it = score_list.begin(), e_it = score_list.end() ;
				it != e_it ; ++it ) {
		std::cout << "Including " << name_from_score_type( *it ) << std::endl;
	}

	for( ScoreTypes::iterator it = score_list.begin(), e_it = score_list.end() ;
				it != e_it ; ++it ) {
		if( std::find( fixed_score_list.begin(), fixed_score_list.end(),  *it ) != fixed_score_list.end() ) {
			// Check the syntax here
			score_list.erase( it );
			std::cout << "Removing " << name_from_score_type( *it ) << std::endl;
		}
	}

	for( ScoreTypes::const_iterator it = score_list.begin(), e_it = score_list.end() ;
				it != e_it ; ++it ) {
		std::cout << "Revised Including " << name_from_score_type( *it ) << std::endl;
	}

	// read list file of pdbs
	utility::vector1< std::string > filenames;

	{
		std::string const listfile( start_file() );
		std::ifstream data( listfile.c_str() );
		std::string line;
		while ( getline( data,line ) ){
			filenames.push_back( line );
		}
		data.close();
	}

	OptEData full_data_set;

	for ( Size n=1; n<= filenames.size(); ++n ) {
		std::string const filename( filenames[n] );

		//std::cout << "Working on file: " << filename << std::endl;

		Pose pose;
		Pose native_pose;

		core::import_pose::centroid_pose_from_pdb( pose, filename , core::import_pose::PDB_file);
		native_pose = pose;

		//std::cout << "read file: " << filename << ' '<< pose.total_residue() << std::endl;

		//std::cout << "SEQUENCE: " << filename << ' ' << pose.sequence() << std::endl;

		//std::cout << "Core optimization version" << std::endl;

		scorefxn( pose );

		//std::cout << "Initial score " << score_orig << std::endl;

		get_opte_data( pose, native_pose, scorefxn, score_list, fixed_score_list, full_data_set );

		//std::cout << "Have data from " << full_data_set.positions() << " positions " << std::endl;
//		utility_exit_with_message( "DONE AFTER FIRST FILE" );
	}

	// Report energies
#ifdef NOTDEF
	for( OptEPositionDataOPs::iterator itr = full_data_set.position_data_begin(),
				e_itr = full_data_set.position_data_end() ; itr != e_itr ; ++itr) {
		std::cout << "Reporting for position " << (*itr)->position() << std::endl;
		RotamerDataOPs &  rot_line( (*itr)->data() );
		for( Size ll = 1, e_ll = rot_line.size() ; ll <= e_ll ; ++ll ) {
			// Report for one rotamer
			std::cout << (*itr)->position() << "\t";
			std::cout << rot_line[ll]->this_aa() << "\t" << rot_line[ll]->rot_number() << "\t";
			for( Size enrgy = 1, e_enrgy = score_list.size() ; enrgy <= e_enrgy ; ++enrgy ) {
			std::cout << (rot_line[ ll ])->[ enrgy ] << "\t";
			}
			std::cout << std::endl;
		}
	}
#endif

	// Do the actual weight optimization
	utility::vector1< Real > component_weights( protocols::optimize_weights::n_optE_data_types, 1.0 );
	OptEMultifunc opt_min( full_data_set, fixed_terms, free_count, score_list, fixed_score_list, component_weights );
	optimization::Multivec start_dofs  = opt_min.get_dofs_from_energy_map( include_terms );

	optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.0000001, true, false );

	optimization::Minimizer minimizer( opt_min, options );
	minimizer.run( start_dofs );

	for ( Size ii = 1 ; ii <= Size( free_count ) ; ++ii ) {
		std::cout << "Weight for " << name_from_score_type( score_list[ ii ] ) << " is " << start_dofs[ ii ] << std::endl;
	}

	for  ( Size ii = free_count+1 ; ii <= Size( free_count+20 ) ; ++ii ) {
		std::cout << "Reference for AA " << ii-free_count << " is " <<  start_dofs[ ii ] << std::endl;
	}

	/// Weight file output:
	std::cout << "METHOD_WEIGHTS ref ";
	for ( Size ii = 1; ii <= Size( free_count+20 ) ; ++ii ) {
		std::cout <<  start_dofs[ ii ] << " ";
	}
	std::cout << std::endl;

	for( Size ii = 1 ; ii <= Size( free_count ) ; ++ii ) {
		std::cout << name_from_score_type( score_list[ ii ] ) << " " << start_dofs[ ii ] << std::endl;
	}


}

//		optimize_weights(
//			pdbs_this_round,
//			native_pdb_names,
//			free_parameters,
//			fixed_parameters,
//			ii_optimized_free_parameters,
//			ii_reference_energies);

void
optimize_weights(
	utility::vector1< std::string > const & filenames,
	utility::vector1< std::string > const & native_filenames,
	core::scoring::EnergyMap const & free_parameters,
	core::scoring::EnergyMap const & fixed_parameters,
	utility::vector1< Real > const & starting_reference_energies,
	core::scoring::EnergyMap & optimized_free_parameters,
	utility::vector1< Real > & optimized_reference_energies
)
{
	using namespace pack::rotamer_set;
	using namespace protocols::optimize_weights;

	core::scoring::EnergyMap include_terms = free_parameters;
	include_terms += fixed_parameters;

	Size include_count( 0 ), fixed_count( 0 ), free_count( 0 );

	ScoreFunctionOP scorefxn = new ScoreFunction;
	ScoreTypes free_score_list;
	for( int i=1 ; i <= n_score_types ; ++i ) {
		if( include_terms[ ScoreType(i) ] != 0.0 ) {
			if( fixed_parameters[ ScoreType(i) ] == 0.0 ) {
				free_score_list.push_back( ScoreType(i) );
				++free_count;
			}
			scorefxn->set_weight( ScoreType(i), include_terms[ ScoreType(i) ] );
			++include_count;
		}
	}

	ScoreTypes fixed_score_list;
	for( int i=1 ; i <= n_score_types ; ++i ) {
		if( fixed_parameters[ ScoreType(i) ] != 0.0 ) {
			fixed_score_list.push_back( ScoreType(i) );
			++fixed_count;
		}
	}

	assert( free_count == include_count - fixed_count );

	OptEData full_data_set;

	for ( Size n=1; n<= filenames.size(); ++n ) {
		std::string const filename( filenames[n] );
		std::string const native_filename( native_filenames[n] );

		//std::cout << "Working on file: " << filename << std::endl;

		core::pose::Pose pose, native_pose;

		core::import_pose::centroid_pose_from_pdb( pose, filename , core::import_pose::PDB_file);
		core::import_pose::centroid_pose_from_pdb( native_pose, native_filename , core::import_pose::PDB_file);

		std::cout << "read file: " << filename << " " << pose.total_residue() << " and native file: " << native_filename << std::endl;

		//std::cout << "SEQUENCE: " << filename << ' ' << pose.sequence() << std::endl;

		//std::cout << "Core optimization version" << std::endl;

		(*scorefxn)( pose );

		//std::cout << "Initial score " << score_orig << std::endl;

		get_opte_data( pose, native_pose, *scorefxn, free_score_list, fixed_score_list, full_data_set );

		//std::cout << "Have data from " << full_data_set.positions() << " positions " << std::endl;
//		utility_exit_with_message( "DONE AFTER FIRST FILE" );
	}

	// Report energies
#ifdef NOTDEF
	for( OptEPositionDataOPs::iterator itr = full_data_set.position_data_begin(),
				e_itr = full_data_set.position_data_end() ; itr != e_itr ; ++itr) {
		std::cout << "Reporting for position " << (*itr)->position() << std::endl;
		RotamerDataOPs & rot_line( (*itr)->data() );
		for( Size ll = 1, e_ll = rot_line.size() ; ll <= e_ll ; ++ll ) {
			// Report for one rotamer
			std::cout << (*itr)->position() << "\t";
			std::cout << rot_line[ll]->this_aa() << "\t" << rot_line[ll]->rot_number() << "\t";
			for( Size enrgy = 1, e_enrgy = free_score_list.size() ; enrgy <= e_enrgy ; ++enrgy ) {
				std::cout << (rot_line[ ll ])->[ enrgy ] << "\t";
			}
			std::cout << std::endl;
		}
	}
#endif

	// Do the actual weight optimization
	utility::vector1< Real > component_weights( protocols::optimize_weights::n_optE_data_types, 1.0 );
	OptEMultifunc opt_min(
		full_data_set, fixed_parameters,
		(int) free_count, free_score_list,
		fixed_score_list, starting_reference_energies,
		component_weights);
	optimization::Multivec start_dofs  = opt_min.get_dofs_from_energy_map( include_terms );

	optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.0000001, true, false );

	optimization::Minimizer minimizer( opt_min, options );
	minimizer.run( start_dofs );

	for ( Size ii = 1 ; ii <= Size( free_count ) ; ++ii ) {
		std::cout << "Weight for " << name_from_score_type( free_score_list[ ii ] ) << " is " << start_dofs[ ii ] << std::endl;
	}

	for  ( Size ii = free_count+1 ; ii <= Size( free_count+20 ) ; ++ii ) {
		std::cout << "Reference for AA " << ii-free_count << " is " <<  start_dofs[ ii ] << std::endl;
	}

	/// Weight file output:
	std::cout << "METHOD_WEIGHTS ref ";
	for ( Size ii = 1; ii <= 20; ++ii ) {
		std::cout <<  start_dofs[ free_count + ii ] << " ";
		optimized_reference_energies[ ii ] = start_dofs[ free_count+ii ]; // save them in non-negated form
	}
	std::cout << std::endl;

	for( Size ii = 1 ; ii <= free_count; ++ii ) {
		std::cout << name_from_score_type( free_score_list[ ii ] ) << " " << start_dofs[ ii ] << std::endl;
	}


	optimized_free_parameters = opt_min.get_energy_map_from_dofs( start_dofs );
	for( Size ii = 1 ; ii <= fixed_score_list.size(); ++ii ) {
		optimized_free_parameters[ fixed_score_list[ ii ]] = 0; // reset fixed parameters
	}
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( optimized_free_parameters[ (ScoreType) ii ] != 0 ) {
			std::cout << "optimized_free_parameters: " << name_from_score_type( (ScoreType) ii ) << " " << optimized_free_parameters[ (ScoreType) ii ] << std::endl;
		}
	}

}
///////////////////////////////////////////////////////////////////////////////

void
iterative_optE()
{
	using namespace core::scoring;
	using namespace core;

	EnergyMap free_parameters;
	EnergyMap fixed_parameters;

	initialize_free_and_fixed( free_parameters, fixed_parameters );

	utility::vector1< std::string > native_pdb_names = get_native_pdb_names();

	utility::vector1< Real > reference_energies( chemical::num_canonical_aas );

	Real last_sequence_recovery_rate( 0.0 );
	utility::vector1< std::string > pdbs_this_round( native_pdb_names );

	for ( Size ii = 1; ii <= 10; ++ii ) {
		bool weights_have_converged( false );

		utility::vector1< std::string > pdbs_next_round;
		setup_pdbnames_next_round( ii, pdbs_next_round, native_pdb_names );

		utility::vector1< Real > ii_reference_energies( reference_energies );

		EnergyMap ii_optimized_free_parameters;

		std::cout << "Round " << ii << std::endl;

		optimize_weights(
			pdbs_this_round,
			native_pdb_names,
			free_parameters,
			fixed_parameters,
			reference_energies,
			ii_optimized_free_parameters,
			ii_reference_energies);

		for ( Size jj = 1; jj <= opte_num_inner_iterations(ii) ; ++jj ) {
			Real jj_scale_factor = opte_weight_mixing_factor( ii, jj );
			std::cout << "Inner Loop round " << jj << " jj_scale_factor: " << jj_scale_factor << std::endl;

			EnergyMap jj_free;
			for ( Size kk = 1; kk <= n_score_types; ++kk ) {
				jj_free[ (ScoreType) kk ] += ( 1.0 - jj_scale_factor ) * free_parameters[ (ScoreType) kk ];
				jj_free[ (ScoreType) kk ] += jj_scale_factor * ii_optimized_free_parameters[ (ScoreType) kk ];
			}

			EnergyMap jj_combined_weights( fixed_parameters );
			jj_combined_weights += jj_free;

			utility::vector1< Real > jj_reference_energies( ii_reference_energies.size() );
			for ( Size kk = 1; kk <= ii_reference_energies.size(); ++kk ) {
				jj_reference_energies[ kk ] += ( 1.0 - jj_scale_factor ) * reference_energies[ kk ];
				jj_reference_energies[ kk ] += jj_scale_factor * ii_reference_energies[ kk ];
			}

			ScoreFunctionOP jj_sfxn = create_score_function_from_weights_and_refEs( jj_combined_weights, jj_reference_energies );
			std::cout << "Measuring sequence recovery with weight set: " << std::endl;
			write_parameters_to_std_out( jj_combined_weights, jj_reference_energies );
			Real jj_seq_recovery = measure_sequence_recovery(
				native_pdb_names,
				pdbs_next_round,
				jj_sfxn );

			if ( accept_new_weight_set( ii, jj, jj_seq_recovery, last_sequence_recovery_rate, jj_scale_factor ) ) {
				last_sequence_recovery_rate = jj_seq_recovery;
				if ( ii != 1 && converged( free_parameters, jj_free, reference_energies, jj_reference_energies ) ){
					weights_have_converged = true;
					break;
				}
				free_parameters = jj_free;
				reference_energies = jj_reference_energies;
				pdbs_this_round = pdbs_next_round;
				break;
			}


		}
		if ( weights_have_converged ) {
			std::cout << "Weights have converged!" << std::endl;
			break;
		}

	}


}

Real
measure_sequence_recovery(
	utility::vector1< std::string > const & native_pdb_names,
	utility::vector1< std::string > const & names_for_output_pdbs,
	core::scoring::ScoreFunctionOP sfxn
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	Size n_residues_total( 0 );
	Size n_residues_recovered( 0 );

	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );

	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( sfxn );

	for ( Size ii = 1; ii <= native_pdb_names.size(); ++ii ) {
		/// read the pdb into a pose
		core::pose::Pose pose;
		if ( option[ in::file::centroid_input ] ) {
			core::import_pose::centroid_pose_from_pdb( pose, native_pdb_names[ ii ] , core::import_pose::PDB_file);
		} else {
			core::import_pose::pose_from_file( pose, native_pdb_names[ ii ] , core::import_pose::PDB_file);
		}

		/// record original sequence
		Size const nres = pose.total_residue();
		n_residues_total += nres;
		utility::vector1< chemical::AA > input_sequence( nres );
		for ( Size jj = 1; jj <= nres; ++jj ) { input_sequence[ jj ] = pose.residue(jj).aa(); }

		/// redesign the pose
		pack_mover->apply( pose );

		/// measure seq recov
		for ( Size jj = 1; jj <= nres; ++jj ) {
			if ( input_sequence[ jj ] == pose.residue(jj).aa()) {
				++n_residues_recovered;
			}
		}

		/// write out new pdb
		pose.dump_pdb( names_for_output_pdbs[ ii ] );

	}
	std::cout << "New Sequence Recovery: " << ( static_cast< Real >  (n_residues_recovered)) / n_residues_total << std::endl;
	return ( static_cast< Real >  (n_residues_recovered)) / n_residues_total;
}


Size
opte_num_inner_iterations(
	Size outer_loop_counter
)
{
	if ( outer_loop_counter == 1 ) return 1;
	else return 5;
}

Real
opte_weight_mixing_factor(
	Size outer_loop_counter,
	Size inner_loop_counter
)
{
	if ( outer_loop_counter == 1 ) {
		return 1.0;
	} else if ( inner_loop_counter <= 5 ) {
		return ( 1.0 / ( inner_loop_counter ) );
	} else {
		return 0.1;
	}
}

core::scoring::ScoreFunctionOP
create_score_function_from_weights_and_refEs
(
	core::scoring::EnergyMap const & emap,
	utility::vector1< Real > const & reference_energies
)
{
	using namespace core::scoring;

	static int n_output_scorefiles( 1 );

	std::stringstream instream;
	instream << n_output_scorefiles;
	++n_output_scorefiles;

	std::string scorefile_name = "optE_scorefile_" + instream.str() + ".wts";

	std::ofstream fout( scorefile_name.c_str() );

	/// Weight file output:
	fout << "METHOD_WEIGHTS ref ";
	for ( Size ii = 1; ii <= reference_energies.size(); ++ii ) {
		fout <<  reference_energies[ ii ] << " ";
	}
	fout << "\n";

	for( Size ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
		if ( emap[ ScoreType( ii ) ] != 0 ) {
			fout << name_from_score_type( ScoreType( ii ) ) << " " << emap[ ScoreType( ii ) ] << "\n";
		}
	}
	fout.close();


	ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( scorefile_name );
	return sfxn;
}

void
initialize_free_and_fixed(
	core::scoring::EnergyMap & free_parameters,
	core::scoring::EnergyMap & fixed_parameters
)
{
	free_parameters[ vdw ] = numeric::random::rg().uniform();
	//free_parameters[ env ] = numeric::random::rg().uniform();
	free_parameters[ pair ] = numeric::random::rg().uniform();
	//free_parameters[ cbeta ] = numeric::random::rg().uniform();
	free_parameters[ rama ] = numeric::random::rg().uniform();
	free_parameters[ p_aa_pp ] = numeric::random::rg().uniform();
	free_parameters[ cenpack ] = numeric::random::rg().uniform();

	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_parameters[ (ScoreType) ii ] != 0 ) {
			std::cout << "Initial free_parameters: " << name_from_score_type( (ScoreType) ii ) << " " << free_parameters[ (ScoreType) ii ] << std::endl;
		}
	}

	fixed_parameters[ env ] = 0.4;
}

bool
converged(
	core::scoring::EnergyMap & free_parameters_prev,
	core::scoring::EnergyMap & free_parameters_curr,
	utility::vector1< Real > const & reference_energies_prev,
	utility::vector1< Real > const & reference_energies_curr
)
{
	using namespace core::scoring;

	assert( reference_energies_prev.size() == reference_energies_curr.size() );

	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( std::abs( free_parameters_prev[ (ScoreType) ii ] - free_parameters_curr[ (ScoreType) ii ]) > 0.001 ) {
			return false;
		}
	}
	for ( Size ii = 1; ii <= reference_energies_prev.size(); ++ii ) {
		if ( std::abs( reference_energies_prev[ ii ] - reference_energies_curr[ ii ]) > 0.001 ) {
			return false;
		}
	}

	std::cout << "Converged: " << std::endl;
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_parameters_prev[ (ScoreType) ii ] != 0 || free_parameters_curr[ (ScoreType) ii ] != 0)
			std::cout << name_from_score_type( ScoreType( ii ) ) << " prev " << free_parameters_prev[ (ScoreType) ii ] << " curr " << free_parameters_curr[ (ScoreType) ii ] << std::endl;
	}
	for ( Size ii = 1; ii <= reference_energies_prev.size(); ++ii ) {
		std::cout <<  ii << " prev " << reference_energies_prev[ ii ] << " curr " << reference_energies_curr[ ii ] << std::endl;
	}


	return true;
}

void
write_parameters_to_std_out(
	core::scoring::EnergyMap & free_parameters,
	utility::vector1< Real > const & reference_energies
)
{
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		if ( free_parameters[ (ScoreType) ii ] != 0 )
			std::cout << name_from_score_type( ScoreType( ii ) ) << " " << free_parameters[ (ScoreType) ii ] << std::endl;
	}
	for ( Size ii = 1; ii <= reference_energies.size(); ++ii ) {
		std::cout << "Reference energy for " <<  ii << " " << reference_energies[ ii ] << std::endl;
	}

}


void
setup_pdbnames_next_round(
	Size const outer_loop_counter,
	utility::vector1< std::string  > & pdbs_next_round,
	utility::vector1< std::string > const & native_pdb_names
)
{
	// assumption, native_pdb_names end in ".pdb"
	pdbs_next_round.resize( native_pdb_names.size() );
	std::stringstream stream;
	stream << outer_loop_counter;
	for ( Size ii = 1; ii <= native_pdb_names.size(); ++ii ) {
		std::string native_substr = native_pdb_names[ ii ].substr( 0, native_pdb_names[ ii ].size() - 4 );
		pdbs_next_round[ ii ] = native_substr + "_" + stream.str() + ".pdb";
	}
}

bool
accept_new_weight_set(
	Size outer_loop_counter,
	Size inner_loop_counter,
	Real latest_seq_recovery_rate,
	Real last_sequence_recovery_rate,
	Real inner_loop_scale_factor
)
{
	if ( outer_loop_counter == 1 ) return true;
	if ( inner_loop_counter == 5 ) return true;
	return (latest_seq_recovery_rate > last_sequence_recovery_rate ) || (inner_loop_scale_factor <= 0.1 );
}

int
main( int argc, char * argv [] )
{
	try {
		//using namespace core;
		devel::init( argc, argv );

		//simple_opte_test();
		iterative_optE();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1
	}
	return 0;
}
