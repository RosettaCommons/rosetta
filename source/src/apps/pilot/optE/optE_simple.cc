// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
//#include <core/options/option.hh>

#include <core/types.hh>


#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>

#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>


#include <protocols/optimize_weights/OptEData.hh>

#include <core/graph/Graph.hh>


#include <core/io/pdb/pose_io.hh>


#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <protocols/optimize_weights/OptEMultifunc.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/util.hh>


#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

//#include <devel/dna/util_public.hh>
//#include <protocols/dna/util.hh>
//#include <protocols/dna/classes.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/RotamerDNAHBondFilter.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

// option key includes

#include <core/options/keys/optE.OptionKeys.gen.hh>
#include <core/options/keys/run.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <ObjexxFCL/format.hh>


//silly using/typedef

using namespace core;
using namespace scoring;
using namespace pack;
using namespace optimization;
using namespace ObjexxFCL::format;

using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
get_opte_data(
	pose::Pose & pose,
	ScoreFunction const & unweighted_scorefxn,
	ScoreTypes & score_list,
	ScoreTypes & fixed_score_vec,
	protocols::optimize_weights::OptEData & opte_data )
{
	using namespace pose;
	using namespace rotamer_set;
	using namespace protocols::optimize_weights;

	protocols::dna::RestrictDesignToProteinDNAInterfaceOP dna_inter_restrictor =
		new protocols::dna::RestrictDesignToProteinDNAInterface();

	task::TaskFactory taskfactory;
	taskfactory.push_back( dna_inter_restrictor );
	task::PackerTaskOP main_task = taskfactory.create_task_and_apply_taskoperations( pose );

	unweighted_scorefxn.setup_for_packing( pose, main_task->repacking_residues(), main_task->designing_residues() );

	vector1< bool > task_mask( pose.total_residue(), false );
	for ( Size pos = 1; pos <= pose.total_residue(); ++pos ) {
		if ( pose.residue_type(pos).is_DNA() ) continue; // skip DNA
//		if ( !main_task->nonconst_residue_task( pos ).being_designed() ) continue;
		if ( !main_task->nonconst_residue_task( pos ).being_packed() ) continue;

		task::PackerTaskOP position_task = task::TaskFactory::create_packer_task( pose );
		position_task->initialize_from_command_line();
		position_task->set_bump_check( false );
		task_mask[ pos ] = true;
		position_task->restrict_to_residues( task_mask );
		task_mask[ pos ] = false;
		// deliver a protein-DNA hbonding filter to the PackerTask for filtering ex rotamers
		// (formerly known as 'rotamer explosion')
		protocols::dna::RotamerDNAHBondFilterOP rot_dna_hb_filter =
			new protocols::dna::RotamerDNAHBondFilter();
		position_task->append_rotamer_operation( rot_dna_hb_filter );

		PNatAAOptEPositionDataOP this_pos_data = new PNatAAOptEPositionData;
		this_pos_data->set_position( pos );
		this_pos_data->set_native_aa( pose.residue( pos ).aa() );
		this_pos_data->set_neighbor_count(
			pose.energies().tenA_neighbor_graph().get_node( pos )->num_neighbors_counting_self() );

		graph::GraphCOP packer_neighbor_graph(
			create_packer_graph( pose, unweighted_scorefxn, position_task ) );

		RotamerSetFactory rsf;
		RotamerSetOP rotset = rsf.create_rotamer_set( pose.residue( pos ) );

		rotset->set_resid( pos );
		rotset->build_rotamers( pose, unweighted_scorefxn, *position_task, packer_neighbor_graph );
		unweighted_scorefxn.prepare_rotamers_for_packing( pose, *rotset );

		// First, need a vector of energy maps
		vector1< EnergyMap > emap_vector( rotset->num_rotamers() );

		// Call the new energy map fn
		rotset->compute_one_body_energy_maps(
			pose, unweighted_scorefxn, *position_task, packer_neighbor_graph, emap_vector );

		for ( Size jj = 1; jj <= rotset->num_rotamers(); ++jj ) {

			EnergyMap & emap_total( emap_vector[jj] );

			// Hacky limit for fa_rep
			if ( emap_total[ fa_rep ] > 20.0 ) emap_total[ fa_rep ] = 20.0;

			vector1< Real > energy_info, fixed_energy_info;

			for ( ScoreTypes::iterator score_type_iter = score_list.begin(),
							end_iter = score_list.end() ; score_type_iter != end_iter; ++score_type_iter ) {
				energy_info.push_back( emap_total[ *score_type_iter ] );
			}

			for ( ScoreTypes::iterator score_type_iter = fixed_score_vec.begin(),
							end_iter = fixed_score_vec.end() ; score_type_iter != end_iter; ++score_type_iter ) {
				fixed_energy_info.push_back( emap_total[ *score_type_iter ] );
			}

			PNatAAOptERotamerDataOP new_rot_line =
				new PNatAAOptERotamerData( (*rotset->rotamer( jj )).aa(), jj, energy_info, fixed_energy_info );

			this_pos_data->add_rotamer_line_data( new_rot_line );
		}

	// Done with rotamers for this position, store this position data object
	opte_data.add_position_data( this_pos_data );

	}
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


	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace optimization;
	using namespace rotamer_set;
	using namespace protocols::optimize_weights;

	ScoreFunction user_scorefxn, unweighted_scorefxn;

	if ( option[ optE::weights ].user() ) {
		// use a supplied weight file to determine score types to be included
		// this should include /all/ relevant energy terms, including those to be fixed
		std::string optE_weights( option[ optE::weights ]() );
		std::cout << "Loaded optE weights from " << optE_weights << std::endl;
		user_scorefxn.initialize_from_file( optE_weights );
		ScoreTypes score_types( user_scorefxn.get_nonzero_weighted_scoretypes() );
		for ( ScoreTypes::const_iterator score_type( score_types.begin() );
		      score_type != score_types.end(); ++score_type ) {
			if ( user_scorefxn.has_zero_weight( *score_type ) ) {
				// useful for pesky weights like fa_intra_rep?
				unweighted_scorefxn.set_weight( *score_type, 0. );
			} else {
				unweighted_scorefxn.set_weight( *score_type, 1. );
			}
		}
	} else  {
		// a common (default) scheme
	  unweighted_scorefxn.set_weight( fa_atr, 1. );
	  unweighted_scorefxn.set_weight( fa_rep, 1. );
	  unweighted_scorefxn.set_weight( fa_sol, 1. );
	  unweighted_scorefxn.set_weight( fa_pair, 1. );
	  unweighted_scorefxn.set_weight( hbond_bb_sc, 1. );
	  unweighted_scorefxn.set_weight( hbond_sc, 1. );
	  unweighted_scorefxn.set_weight( fa_dun, 1. );
	  unweighted_scorefxn.set_weight( p_aa_pp, 1. );
	}

	// Make two EnergyMaps - one is the set of weights which are included in the optimization, the other indicates which should be fixed at a certain non-zero value
	EnergyMap include_terms( unweighted_scorefxn.weights() ), fixed_terms;

	if ( option[ optE::fix ].user() ) {
		// command-line fixing of weights (that also appear in the command-line-supplied weight file)
		vector1< std::string > fixed_types( option[ optE::fix ]() );
		for ( vector1< std::string >::const_iterator fixed_type( fixed_types.begin() );
		      fixed_type != fixed_types.end(); ++fixed_type ) {
			ScoreType fixed_score_type( score_type_from_name( *fixed_type ) );
			if ( !user_scorefxn.has_nonzero_weight( fixed_score_type ) ) {
				std::cout << "score type to fix: \"" << fixed_score_type
				          << "\" not found in user-defined scorefxn!" << std::endl;
			}
			fixed_terms.set( fixed_score_type, user_scorefxn.get_weight( fixed_score_type ) );
		}
	} else {
		// a common (default) scheme
		fixed_terms.set( fa_atr, 0.8 );
	}

	// Loop through the energy maps and figure out how many terms are included, and how many total are free.

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
			std::cout << "Removing " << name_from_score_type( *it ) << std::endl;
			score_list.erase( it );
		}
	}

	for( ScoreTypes::const_iterator it = score_list.begin(), e_it = score_list.end() ;
				it != e_it ; ++it ) {
		std::cout << "Revised Including " << name_from_score_type( *it ) << std::endl;
	}

	OptEData full_data_set( fixed_score_list, score_list );

	if ( option[ optE::data_in ].user() ) {
		std::string data_file( option[ optE::data_in ]() );
		std::cout << "Reading optE data from binary file " << data_file << std::endl;
		//full_data_set.read_from_file( data_file );
		full_data_set.read_from_binary_file( data_file );
	} else {

		vector1< std::string > pdb_files( start_files() );
		for ( vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
		      pdb_file != pdb_files.end(); ++pdb_file ) {

			std::cout << "Working on file: " << *pdb_file << std::endl;

			Pose pose;
			core::import_pose::pose_from_pdb( pose, *pdb_file );

			std::cout << "read file: " << *pdb_file << ' '<< pose.total_residue() << std::endl;
			std::cout << "SEQUENCE: " << *pdb_file << ' ' << pose.sequence() << std::endl;

			std::cout << "Core optimization version" << std::endl;
			Energy score_orig = unweighted_scorefxn( pose );
			std::cout << "Initial score " << score_orig << std::endl;

			get_opte_data( pose, unweighted_scorefxn, score_list, fixed_score_list, full_data_set );

			std::cout << "Have data from " << full_data_set.num_rotamers() << " rotamers at "
			          << full_data_set.num_positions() << " positions " << std::endl;
	//		utility_exit_with_message( "DONE AFTER FIRST FILE" );
		}
	}

#ifdef NOTDEF
	// Report energies
	for( OptEPositionDataOPs::iterator itr = full_data_set.position_data_begin(),
				e_itr = full_data_set.position_data_end() ; itr != e_itr ; ++itr) {
		std::cout << "Reporting for position " << (*itr)->position() << std::endl;
		OptERotamerDataOPs &  rot_line(  (*itr)->data() );
		for( Size ll = 1, e_ll = rot_line.size() ; ll <= e_ll ; ++ll ) {
			// Report for one rotamer
			std::cout << (*itr)->position() << "\t";
			std::cout << rot_line[ll].this_aa() << "\t" << rot_line[ll].rot_number() << "\t";
			for( Size enrgy = 1, e_enrgy = score_list.size() ; enrgy <= e_enrgy ; ++enrgy ) {
			std::cout << (rot_line[ ll ])[ enrgy ] << "\t";
			}
			std::cout << std::endl;
		}
	}
#endif

	if ( option[ optE::data_out ].user() ) {
		//full_data_set.write_to_file( option[ optE::data_out ]() );
		full_data_set.write_to_binary_file( option[ optE::data_out ]() );
	}

	// Do the actual weight optimization
	utility::vector1< Real > component_weights( protocols::optimize_weights::n_optE_data_types, 1.0 );
	OptEMultifunc opt_min( full_data_set, fixed_terms, free_count, score_list, fixed_score_list, component_weights );

	optimization::Multivec start_dofs = opt_min.get_dofs_from_energy_map( include_terms );

	for ( optimization::Multivec::const_iterator it( start_dofs.begin() ); it != start_dofs.end();
	      ++it ) {
		std::cout << F(3,2,*it) << " ";
	}
	std::cout << std::endl;

	if ( option[ optE::weights ].user() ) {
		//ja optional: use initial weights from the weights file
		std::string optE_weights( option[ optE::weights ]() );

		// get starting reference weights and inform the OptEMultifunc
		vector1< Real > ref_energies =
			user_scorefxn.energy_method_options().method_weights( score_type_from_name( "ref" ) );
		opt_min.set_starting_reference_energies( ref_energies );

		start_dofs = opt_min.get_dofs_from_energy_map( user_scorefxn.weights() );

		std::cout << "Using starting weights from " << optE_weights << ":" << std::endl;
		for ( optimization::Multivec::const_iterator it( start_dofs.begin() ); it != start_dofs.end();
		      ++it ) {
			std::cout << F(3,2,*it) << " ";
		}
		std::cout << std::endl;
	}

	std::string min_type("dfpmin_armijo_nonmonotone");
	if ( option[ run::min_type ].user() ) min_type = option[ run::min_type ]();
	optimization::MinimizerOptions options( min_type, 0.0000001, true, true );

	optimization::Minimizer minimizer( opt_min, options );
	minimizer.run( start_dofs );

	// output a weight file
	std::string name( "optE_fit.wts" );
	if ( option[ optE::data_in ].user() ) name = option[ optE::data_in ]() + ".wts";
	std::ofstream weight_file( name.c_str() );

	weight_file << "# weights trained by optE" << std::endl;

	weight_file << "# fixed params:";
	for ( ScoreTypes::const_iterator st( fixed_score_list.begin() );
	      st != fixed_score_list.end(); ++st ) {
		weight_file << " " << *st;
	}
	weight_file << std::endl;

	weight_file << "# free params:";
	for ( ScoreTypes::const_iterator st( score_list.begin() );
	      st != score_list.end(); ++st ) {
		weight_file << " " << *st;
	}
	weight_file << std::endl;

	for( int i(1); i <= n_score_types; ++i ) {
		if( fixed_terms[ ScoreType(i) ] != 0.0 ) {
			weight_file << ScoreType(i) << " " << fixed_terms[ ScoreType(i) ] << std::endl;
		}
	}
	Real hbond_bb_sc_weight(0.0);
	for ( Size ii = 1 ; ii <= Size( free_count ) ; ++ii ) {
		std::string type_name( name_from_score_type( score_list[ ii ] ) );
		weight_file << type_name << " " << start_dofs[ ii ] << std::endl;
		if ( type_name == "hbond_bb_sc" ) hbond_bb_sc_weight = start_dofs[ii];
	}
	// hack the hbond_bb weights to be equal to hbond_bb_sc
	if ( hbond_bb_sc_weight != 0.0 ) {
		weight_file << hbond_lr_bb << " " << hbond_bb_sc_weight << std::endl;
		weight_file << hbond_sr_bb << " " << hbond_bb_sc_weight << std::endl;
	}
	weight_file << "ref 1" << std::endl;
	weight_file << "METHOD_WEIGHTS ref";
	for( Size ii = free_count+1 ; ii <= Size( free_count+20 ) ; ++ii ) {
		weight_file << " " << start_dofs[ ii ];
	}
	weight_file << std::endl;
	weight_file.close();
	std::cout << "wrote weights to file " << name << std::endl;

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try{
	//using namespace core;
	devel::init( argc, argv );

	simple_opte_test();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
