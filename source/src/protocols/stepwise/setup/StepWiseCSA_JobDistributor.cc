// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/setup/StepWiseCSA_JobDistributor.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/setup/StepWiseCSA_JobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/EnergyNames.hh>
#include <core/io/silent/util.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/util.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/basic_sys_util.hh>

#include <basic/Tracer.hh>
#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.stepwise.setup.StepWiseCSA_JobDistributor" );

using namespace core;
using namespace core::io::silent;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Initial attempt at population-based optimization for stepwise monte carlo.
//
//  CSA = conformational space annealing. Not attempting to directly use that Scheraga/Lee method,
//         just recalling from memory what M. Tyka and I coded up a while back for Rosetta++.
//
//  nstruct       = # number of updates per number of structures in bank.
//  cycles        = # monte carlo cycles per update (specify as -cycles)
//  csa_bank_size = total # of cycles to carry out by all nodes over all calculation.
//
// So the total compute = total cycles = csa_bank_size * nstruct * cycles.
//
// Instead of a master/worker setup, each node knows the name of the silent file with the
//     "bank" of models and updates it after doing some monte carlo steps. Each node sets a lock
//     by creating a ".lock" file when it needs to read it in and update the bank.
//     Admittedly, kind of janky to use file server for communication, but hopefully will work,
//     and later can fixup for MPI, threading, or whatever.
//
// Each model in a bank has a "cycles" column (and its name is S_N, where N = cycles). This
//     is the total # cycles in the CSA calculation over all nodes completed at the time
//     the model is saved to disk.
//
// Currently decisions to replace models with 'nearby' models are based on RMSD -- note trick
//     below where one pose is "filled" in to be a complete pose. This may take some extra computation,
//     and if the nodes get jammed up waiting for RMSDs to be calculate, this could be optimized (e.g.,
//     by precomputing RMSDs ahead of time, and only locking the file at the last moment.)
//
// Could almost certainly define a "recombination" protocol, and make this into a reasonable
//     genetic algorithm.
//
//      -- rhiju, 2014
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace setup {

	//Constructor
 	StepWiseCSA_JobDistributor::StepWiseCSA_JobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
																													std::string const silent_file,
																													core::Size const nstruct,
																													core::Size const csa_bank_size,
																													core::Real const csa_rmsd,
																													bool const output_round_silent_files ):
		StepWiseJobDistributor( stepwise_monte_carlo, silent_file, nstruct ),
		lock_file( silent_file + ".lock" ),
		csa_bank_size_( csa_bank_size ),
		total_updates_( csa_bank_size * nstruct ),
		csa_rmsd_( csa_rmsd ),
		output_round_silent_files_( output_round_silent_files ),
		total_updates_so_far_( 0 )
	{}

	//Destructor
	StepWiseCSA_JobDistributor::~StepWiseCSA_JobDistributor()
	{}

	void
	StepWiseCSA_JobDistributor::initialize( core::pose::Pose const & pose ) {
		start_pose_ = pose.clone();
		total_updates_so_far_ = 0;
	}

	////////////////////////////////////////////////////////////
	bool
	StepWiseCSA_JobDistributor::has_another_job() {
		return ( total_updates_so_far_ < total_updates_ );
	}

	////////////////////////////////////////////////////////////
	Size
	StepWiseCSA_JobDistributor::get_updates( core::io::silent::SilentStructCOP s ) const {
		runtime_assert( s->has_energy( "updates" ) );
		return static_cast< Size >( s->get_energy( "updates" ) );
	}

	void
	StepWiseCSA_JobDistributor::set_updates( core::io::silent::SilentStructOP s, Size const updates ) const {
		s->add_string_value( "updates", utility::to_string( updates ) );
		s->set_decoy_tag( "S_" +  ObjexxFCL::lead_zero_string_of( updates, 6 ) ); // redundant with updates column.
	}

	////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::apply( core::pose::Pose & pose ) {

		///////////////////////
		// "check out" a model
		///////////////////////
		if ( sfd_ == 0 ) read_in_silent_file();
		if ( !has_another_job() ) return;

		////////////////////
		// run some cycles
		////////////////////
		if ( sfd_->structure_list().size() < nstruct_ ) {
			// If bank is not yet full, start from scratch.
			pose = *start_pose_;
			stepwise_monte_carlo_->set_model_tag( "NEW" );
		} else {
			// Start from a model in the bank.
			SilentStructOP s = numeric::random::rg().random_element( sfd_->structure_list() );
			TR << TR.Cyan << "Starting from model in bank " << s->decoy_tag() <<  TR.Reset << std::endl;
			s->fill_pose( pose );
			stepwise_monte_carlo_->set_model_tag( s->decoy_tag() );
		}
 		stepwise_monte_carlo_->apply( pose );

		/////////////////////
		// now update bank
		/////////////////////
		put_lock_on_silent_file();
		read_in_silent_file();
		update_bank( pose );
		write_out_silent_file();
		free_lock_on_silent_file();

		if ( output_round_silent_files_ ) write_out_round_silent_file();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::update_bank( core::pose::Pose & pose ) {

		runtime_assert( sfd_ != 0 );

		SilentStructOP s; // this eventually most hold the silent struct that is *updated*
		// creates a mock version of the pose with all residues instantiated. used anyway for rms_fill calculation,
		// and below used to check for closeness with poses in bank.
		core::pose::PoseOP full_model_pose;

		if ( sfd_->size() >= csa_bank_size_ ) {
			// not in initial bank-filling phase -- need to make a choice about what to keep in bank
			runtime_assert( sfd_->size() == csa_bank_size_ ); // sanity check.

			// need to decide whether to insert into bank or not.
			sfd_->order_by_energy();
			SilentStructOP s_worst = sfd_->structure_list()[ sfd_->size() ]; // note that sfd readin always orders by energy.
			Real worst_score = s_worst->get_energy( "score" );
			Real score = pose.energies().total_energy();

			// if worse in energy than everything else, don't save. but record # updates in lowest scoring pose.
			if ( score > worst_score ) {
				s = sfd_->structure_list()[ 1 ];
			} else {
				// else if better in energy ... there's a chance this will get inserted.
				utility::vector1< SilentStructOP > const & struct_list = sfd_->structure_list();
				full_model_pose = monte_carlo::build_full_model( pose ); // makes a full model with all residues.

				//  If within X RMSD of a model, replace it.
				Size kick_out_idx( struct_list.size() );
				for ( Size n = 1; n <= struct_list.size(); n++ ){
				 	SilentStructOP s_test = struct_list[ n ];
				 	Pose pose_test;
				 	s_test->fill_pose( pose_test );
				 	if ( check_for_closeness( pose_test, *full_model_pose ) ){
						Real score_test = s_test->get_energy( "score" );
						TR << TR.Magenta << "Found a pose that was close. " << score_test << " vs. " << score << TR.Reset << std::endl;
						if ( score_test > score ) {
							kick_out_idx = n; // replace the pose in the bank with this one.
						} else {
							kick_out_idx = 0;
							s = s_test; // keep the lower energy pose that is already in the bank
						}
						break;
				 	}
				}

				if ( kick_out_idx > 0 ){
					runtime_assert( s == 0 ); // s will be filled below.
					SilentFileDataOP sfd_new( new SilentFileData );
					for ( Size n = 1; n <= struct_list.size(); n++ ) if ( n != kick_out_idx ) sfd_new->add_structure( struct_list[ n ] );
					sfd_ = sfd_new;
				}
			}
		}

		if ( s == 0 ) {
			s = monte_carlo::prepare_silent_struct( "S_0", pose, get_native_pose(),
									superimpose_over_all_, true /*do_rms_fill_calculation*/, full_model_pose );
			sfd_->add_structure( s );
		}

		// allows for checks on book-keeping
		runtime_assert( s != 0 );
		total_updates_so_far_ += 1;
		set_updates( s, total_updates_so_far_ );
		TR << TR.Cyan << "Outputting silent structure: " << s->decoy_tag() << TR.Reset << std::endl;
		sfd_->order_by_energy();

	}

	/////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::read_in_silent_file(){

		sfd_ = SilentFileDataOP( new SilentFileData );
		total_updates_so_far_ = 0;
		if ( !utility::file::file_exists( silent_file_ ) ) return;

		sfd_->read_file( silent_file_ );
		sfd_->order_by_energy();
		utility::vector1< SilentStructOP > const & struct_list = sfd_->structure_list();
		for ( Size n = 1; n <= struct_list.size(); n++ ){
			SilentStructOP s = struct_list[ n ];
			total_updates_so_far_ = std::max( total_updates_so_far_, get_updates( s ) );
		}

	}

	/////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::write_out_silent_file( std::string const silent_file_in ){
		std::string const silent_file = silent_file_in.size() == 0 ?  silent_file_ : silent_file_in;
		runtime_assert( sfd_ != 0 );
		runtime_assert( sfd_->structure_list().size() > 0 );
		runtime_assert( sfd_->structure_list().size() <= csa_bank_size_ );
		stepwise::modeler::remove_silent_file_if_it_exists( silent_file );
		runtime_assert( !utility::file::file_exists( silent_file ) );
		sfd_->clear_shared_silent_data(); // kind of bananas, but otherwise getting header printed out twice and lots of issues.
		sfd_->write_all( silent_file );

	}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::write_out_round_silent_file(){
		// keep track of stages along the way -- every time we do a multiple of csa_bank_size is a 'round'.
		if ( ( total_updates_so_far_ % csa_bank_size_ ) == 0 ){
			Size const nrounds = total_updates_so_far_ / csa_bank_size_;
			std::string const round_silent_file = utility::replace_in( silent_file_, ".out",
          ".round" + ObjexxFCL::lead_zero_string_of( nrounds, 3 ) + ".out" );
			write_out_silent_file( round_silent_file );
			TR << "Done with output: " << round_silent_file << std::endl;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::put_lock_on_silent_file() {
		Size ntries( 0 );
		Size const MAX_TRIES( 1000 );
		while ( utility::file::file_exists( lock_file ) && ++ntries <= MAX_TRIES ){
			TR << lock_file << " exists. Will try again in one second." << std::endl;
			utility::rand_sleep();
		}
		runtime_assert( ntries <= MAX_TRIES );
		std::ofstream ostream( lock_file.c_str() );
		ostream << "locking " << silent_file_ << " by creating: " << lock_file << std::endl;
		ostream.close();
	}

	////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseCSA_JobDistributor::free_lock_on_silent_file() {
		runtime_assert( utility::file::file_exists( lock_file ) );
		utility::file::file_delete( lock_file );
	}


	////////////////////////////////////////////////////////////////////////////////////////
	// The reason to use full_model_pose, which has all residues filled in with A-form, is
	//  that we don't need to do any crazy book-keeping to match the various sister poses in one
	//  model to sister poses in the other model.
	//
	// Also we create full_model_pose anyway to calculate rms_fill.
	//
	bool
	StepWiseCSA_JobDistributor::check_for_closeness( pose::Pose & pose_test,
																									 pose::Pose const & full_model_pose ) const {
		Real const rms = modeler::align::superimpose_with_stepwise_aligner( pose_test, full_model_pose, superimpose_over_all_ );
		TR << TR.Cyan << "Calculated RMS to model: " << tag_from_pose( pose_test ) << " to be: " << rms << " and compared to " << csa_rmsd_ << TR.Reset << std::endl;
		return ( rms < csa_rmsd_ );
	}


} //setup
} //stepwise
} //protocols
