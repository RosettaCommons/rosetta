// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Emperor.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/mpi_refinement/MPI_Refinement.hh>
#include <protocols/mpi_refinement/MPI_Refine_Emperor.hh>
#include <protocols/wum/WorkUnitBase.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <ObjexxFCL/format.hh>

/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

#ifndef _WIN32 // REQUIRED FOR WINDOWS

#include <utility/vector1.hh>

#endif

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace mpi_refinement {

using namespace protocols::wum;

static basic::Tracer TR("MPI.LHR.E");

void
MPI_Refine_Emperor::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_max_lib_size( option[ OptionKeys::lh::max_emperor_lib_size]() );
	max_emperor_lib_round_ = option[ OptionKeys::lh::max_emperor_lib_round ]();
	dump_rounds_ = option[ OptionKeys::lh::write_all_fa_structs ]();

	// report frequency
	n_accept_cummul_ = 0;
	n_change_emperor_report_ = 5; // Report if library added/replaced more than this number
	n_addcall_ = 0;
	n_addcall_emperor_report_ = 1; // Report if add_structures called more than this times
	n_dump_ = 0;
	process_termination_ = false; // super-switch when terminating everything
	termination_broadcasted_ = false;
}


void
MPI_Refine_Emperor::init(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Are we resuming an old job ?
	if ( mpi_resume() != "" ) {
		TR << "Resuming job from IDENT:  " <<  mpi_resume() << std::endl;
		load_state( mpi_resume() );
	};

	if ( option[ OptionKeys::lh::mpi_read_structure_for_emperor ]() ) {
		load_structures_from_cmdline_into_library( library_ref() );
	}

	// Emperors have different library rules!
	set_mpi_feedback( "add_n_replace" );

	print_library( library_central(), "EmpLIB Init ");
}

void
MPI_Refine_Emperor::go()
{
	// initialize master (this is a virtual functino call and this function is overloaded by the children of this class)
	TR << "Init Emperor on mpi_rank " << mpi_rank() << std::endl;
	init();

	TRDEBUG << "Emperor Node: Waiting for data ..." << std::endl;
	while ( true ) {
		// process any incoming messages such as incoming
		TRDEBUG << "Emperor: processing msgs.." << std::endl;
		process_incoming_msgs();

		TRDEBUG << "Emperor: process incoming: " << inbound().size() << std::endl;
		process_inbound_wus();  // lets borrow a master's routine

		// Send to all masters and breakdown itself immediately
		TRDEBUG << "Emperor: processing termination.." << std::endl;
		bool terminate = process_termination();
		if ( terminate ) break;

		TRDEBUG << "Emperor: process outbound" << std::endl;
		process_outbound_wus();// lets borrow a master's routine

		TRDEBUG << "stucked here?" << std::endl;
		// ok, we've done all our work, now wait until we hear from our slaves/masters
		process_incoming_msgs( true );
		TRDEBUG << "not stucked here!" << std::endl;

		print_stats_auto();

	}

	TR << "Process terminating!" << std::endl;
	TR << "Report Total time: " << std::endl;
	report_time();
	// Report final structures at
	dump_structures( library_central(), false, "final" );
}

void
MPI_Refine_Emperor::process_inbound_wus(){

	if ( inbound().size() > 0 ) {
		TRDEBUG << "Processing inbound WUs on emperor .." << std::endl;
	}

	while ( inbound().size() > 0 )
			{
		TRDEBUG << "inbound? " << inbound().next()->get_wu_type() << std::endl;
		WorkUnitBaseOP  next_wu =  inbound().pop_next();
		runtime_assert( next_wu );
		WorkUnit_SilentStructStoreOP structure_wu =
			utility::pointer::dynamic_pointer_cast<  WorkUnit_SilentStructStore > ( next_wu );

		if ( structure_wu.get() == nullptr ) {
			TR << "Cannot save structural data for WU: " << std::endl;
			next_wu->print( TR );
			continue;
		}

		SilentStructStore &decoys = structure_wu->decoys();
		decoys.all_sort_silent_scores();

		// Register master rank so as to moniter procedures...
		core::Size master_rank( structure_wu->last_received_from() );
		if ( masters_done_.find( master_rank ) == masters_done_.end() ) {
			masters_done_[master_rank] = false;
		}

		// Master -> emperor ( just add )
		if ( structure_wu->get_wu_type() == "resultstore" ) {
			TR << "Emperor: received structures: " << decoys.size() << std::endl;
			add_structures_to_library( decoys, "NSGAII" );

			// send back anything
			//WorkUnit_WaitOP wait_wu = new WorkUnit_Wait();
			//wait_wu->set_wu_type( "waitwu" );
			//send_MPI_workunit( wait_wu, structure_wu->last_received_from() );

		} else if ( structure_wu->get_wu_type() == "resultfeedback" ) {
			// Master -> Emperor ( add and request for new structure )
			core::Size const nreceived( decoys.size() );
			TR << "Emperor: received structures: " << nreceived << std::endl;

			// Add via genetic algorithm !
			add_structures_to_library( decoys, "NSGAII" );

			// use extra_data_3 as nsend
			core::Size nsend_requested = structure_wu->extra_data_3();
			core::Size nsend = ( nsend_requested < library_central().size())?
				nsend_requested : library_central().size();

			//std::string const send_option( structure_wu->get_options() );
			std::string const objfunction( structure_wu->get_options() );
			std::string pick_strategy( "random" );
			if ( structure_wu->extra_data_2() == 1 ) {
				pick_strategy = "weighted";
			} else if ( structure_wu->extra_data_2() == 2 ) {
				pick_strategy = "sort";
			}

			// Send back
			TR << "Sending " << nsend << " new structure to master " << structure_wu->last_received_from();
			TR << ", based on objfunction " << objfunction << std::endl;

			if ( pick_strategy.compare( "sort" ) == 0 && nsend_requested <= library_central().size() ) {
				send_sortedpick_library_structs( structure_wu->last_received_from(),
					nsend, objfunction, false );
			} else if ( pick_strategy.compare( "weighted" ) == 0 && nsend_requested <= library_central().size() ) {
				send_sortedpick_library_structs( structure_wu->last_received_from(),
					nsend, objfunction, true );
			} else {
				send_random_library_structs( structure_wu->last_received_from(), nsend );
			}
			TR.Debug << "Sending done." << std::endl;

		} else if ( structure_wu->get_wu_type() == "getnewstruct" ) {
			// Note: Below is not going to be used in Genetic Algorithm.
			TR << "Emperor: received signal for new structure " << decoys.size() << std::endl;
			// dump structures
			if ( decoys.size() > 0 ) {
				add_structures_to_library( decoys );

				// Always send back a structure so the master can continue its life.
				TR << "Sending a new random structure to master:" << structure_wu->last_received_from() << std::endl;

				// send back a random ssid. crude way to prevent interference from old returning structures on the master.
				// really crude. sorry, this could be better i know.
				send_random_library_struct( structure_wu->last_received_from(),
					(core::Size) numeric::random::random_range(0,9999) );
				TR << "Done." << std::endl;
			}

			// Termination from master:
			// collect signal from all masters in order to terminate whole processes
		} else if ( structure_wu->get_wu_type() == "terminate" ) {

			core::Size imaster = structure_wu->last_received_from();
			TR << "Received termination signal from Master no. " << imaster << std::endl;

			// Send back anything to the waiting Master
			WorkUnit_WaitOP wait_wu( new WorkUnit_Wait() );
			wait_wu->set_wu_type( "waitwu" );
			send_MPI_workunit( wait_wu, structure_wu->last_received_from() );

			process_termination_ = true;

		} else if ( structure_wu->get_wu_type() == "terminated" ) {

			core::Size imaster = structure_wu->last_received_from();
			TR << "Confirmed termination on Master no. " << imaster << std::endl;
			masters_done_[imaster] = true;

		} else {
			TR.Error << "ERROR: Unknown workunit received. " << std::endl;
		}

	}

	save_state_auto();
	print_stats();
}

void
MPI_Refine_Emperor::process_outbound_wus(){
}

bool
MPI_Refine_Emperor::process_termination(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::Size nmasters_tot( option[ OptionKeys::wum::n_masters ]() );

	if ( process_termination_ ) {
		// Erase all the queues first
		for ( auto iter = outbound().begin(); iter != outbound().end(); ) {
			iter->reset();
			iter = outbound().erase( iter );
		}

		if ( !termination_broadcasted_ ) {
			WorkUnit_SilentStructStoreOP terminate_wu( new WorkUnit_SilentStructStore() );
			terminate_wu->set_wu_type( "terminate" );

			// 0 is emperor itself
			// masters
			for ( int i = 1; i <= (int)(nmasters_tot); ++i ) {
				//if( !masters_done_[i] ){
				TR << "Sent termination signal to master " << i << std::endl;
				send_MPI_workunit( terminate_wu, i ); // The 0 is the MPI_RANK of the emperror
			}
			termination_broadcasted_ = true;
		}
	}

	core::Size nmaster_done( 0 );
	for ( std::map< core::Size, bool >::const_iterator it = masters_done_.begin();
			it != masters_done_.end(); ++it ) {
		if ( it->second ) nmaster_done++;
	}

	//TR << "Finished " << nmaster_done << " so far out of " << nmasters_tot << std::endl;

	// Alive until all the masters are terminated
	if ( nmaster_done >= nmasters_tot ) {
		return true;
	} else {
		return false;
	}
}

// override the MPI_Refinement function
bool
MPI_Refine_Emperor::add_structures_to_library( SilentStructStore &new_structs,
	std::string add_algorithm ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool result = false;

	SilentStructStore accepted_structs;

	// number of this function calls
	n_addcall_++;

	if ( add_algorithm.compare("NSGAII") == 0 ) {
		// Library_central will be updated, removed structs will be returned in new_structs
		//core::Real const simlimit = option[ OptionKeys::lh::rms_limit ]();
		result = fobj_->update_library_NSGAII( library_central(), new_structs, max_lib_size(), true );

	} else {
		for ( SilentStructStore::const_iterator it = new_structs.begin();
				it != new_structs.end(); ++it ) {
			runtime_assert( *it );
			core::io::silent::SilentStructOP pss = *it;

			// Filter for max_emperor_lib_round_
			if ( max_emperor_lib_round_ > 0 ) {
				core::Size structure_round = (core::Size) pss->get_energy("round");
				if ( structure_round > max_emperor_lib_round_ ) continue;
			}

			// add the structure if it passes energy and rms filters evaluated further down there
			runtime_assert( pss );
			pss->add_energy( "emperor_count",  pss->get_energy( "emperor_count" ) );
			bool local_result = add_structure_to_library( pss, add_algorithm );
			result |= local_result;
			if ( local_result ) {
				accepted_structs.add( *pss );
				n_accept_cummul_ ++;
			}
		}
		if ( result ) limit_library();
	}

	library_central().sort_by( "score" );
	library_central().sort_by( "frontier" );

	std::stringstream prefix("");
	prefix << "iter" << n_addcall_;
	retag_library( library_central(), prefix.str()  );
	library_central().all_add_energy( "iter", n_addcall_ );

	if ( n_accept_cummul_ >= n_change_emperor_report_ ) {
		TR << "=========================================================" << std::endl;
		TR << "Reporting Emperor library based on library change at time:" << std::endl;
		TR << "=========================================================" << std::endl;
		print_library( library_central(), "EmpLIB " );
		TR << "---------------------------------------------------------" << std::endl;
		TR << "removed structures:" << std::endl;
		print_library( new_structs, "EmpDEL " );
		TR << "=========================================================" << std::endl;

		n_accept_cummul_ = 0;

		print_summary( "SUMME "  );

	} else if ( n_addcall_%n_addcall_emperor_report_ == 0  ) {
		TR << "=========================================================" << std::endl;
		TR << "Reporting Emperor library based on library addition call." << std::endl;
		TR << "=========================================================" << std::endl;
		print_library( library_central(), "EmpLIB " );
		TR << "---------------------------------------------------------" << std::endl;
		TR << "removed structures:" << std::endl;
		print_library( new_structs, "EmpDEL " );
		TR << "=========================================================" << std::endl;

		print_summary( "SUMME "  );
	}

	// 1. Always dump *everything* returned to emperor!
	if ( dump_rounds_ ) {
		n_dump_++;
		dump_structures( library_central(), false, "round"+string_of( n_dump_ )+"_" );
	}

	return result;
}


} // namespace mpi_refinement
} // namespace protocols


