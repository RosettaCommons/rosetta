// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Emperor.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/loophash/MPI_LoopHashRefine.hh>
#include <protocols/loophash/MPI_LoopHashRefine_Emperor.hh>
// AUTO-REMOVED #include <protocols/loophash/WorkUnit_LoopHash.hh>
#include <protocols/wum/WorkUnitBase.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <ObjexxFCL/format.hh>
/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

#ifndef _WIN32 // REQUIRED FOR WINDOWS
// AUTO-REMOVED #include <unistd.h>

#include <utility/vector1.hh>

#endif

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loophash {

using namespace protocols::wum;

static basic::Tracer TR("MPI.LHR.Emperor");

static numeric::random::RandomGenerator RG(1248321);  // <- Magic number, do not change it (and dont try and use it anywhere else)


void
MPI_LoopHashRefine_Emperor::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_max_lib_size( option[ OptionKeys::lh::max_emperor_lib_size]() );
	max_emperor_lib_round_ = option[ OptionKeys::lh::max_emperor_lib_round ]();
}


void
MPI_LoopHashRefine_Emperor::init(){
	// Are we resuming an old job ?
	if( mpi_resume() != "" ){
		TR << "Resuming job from IDENT:  " <<  mpi_resume() << std::endl;
		load_state( mpi_resume() );
	};

	// Emperors have different library rules!
	set_mpi_feedback( "add_n_replace" );

	TR << "STARTLIB: " << std::endl;
	print_library();
}

void
MPI_LoopHashRefine_Emperor::go()
{
	// initialize master (this is a virtual functino call and this function is overloaded by the children of this class)
	TR << "Init Master: " << mpi_rank() << std::endl;
	init();

	TR << "Emperor Node: Waiting for data ..." << std::endl;
	while(true){
		// process any incoming messages such as incoming
		TRDEBUG << "Emperor: processing msgs.." << std::endl;
		process_incoming_msgs();

		TRDEBUG << "Emperor: process incoming" << std::endl;
		process_inbound_wus();  // lets borrow a master's routine

		TRDEBUG << "Emperor: process outbound" << std::endl;
		process_outbound_wus();// lets borrow a master's routine

		// ok, we've done all our work, now wait until we hear from our slaves/masters
		process_incoming_msgs( true );

		print_stats_auto();
	}
}


void
MPI_LoopHashRefine_Emperor::process_inbound_wus(){
	if( inbound().size() > 0 ){
		TR << "Processing inbound WUs on emperor .." << std::endl;
	}

	while( inbound().size() > 0 )
	{
		WorkUnitBaseOP  next_wu =  inbound().pop_next();
		runtime_assert( next_wu );
		WorkUnit_SilentStructStoreOP structure_wu = dynamic_cast<  WorkUnit_SilentStructStore * > ( next_wu() );

		if ( structure_wu.get() == NULL ){
			TR << "Cannot save structural data for WU: " << std::endl;
			next_wu->print( TR );
			continue;
		}

		SilentStructStore &decoys = structure_wu->decoys();
		decoys.all_sort_silent_scores();
		if ( structure_wu->get_wu_type() == "resultpack" ){
			TR << "Emperor: receivd structures: " << decoys.size() << std::endl;
			// dump structures
			add_structures_to_library( decoys );
		} else
		if ( structure_wu->get_wu_type() == "getnewstruct" ){
			TR << "Emperor: received expired structures: " << decoys.size() << std::endl;
			// dump structures
			if( decoys.size() > 0 ){
				add_structures_to_library( decoys );
				// Always send back a structure so the master can continue its life.
				TR << "Sending a new random structure to master:" << structure_wu->last_received_from() << std::endl;
				// send back a random ssid. crude way to prevent interference from old returning structures on the master. really crude. sorry, this could be better i know.
				send_random_library_struct( structure_wu->last_received_from(), (core::Size) numeric::random::random_range(0,9999) );
				TR << "Done." << std::endl;
			}
		} else {
			TR.Error << "ERROR: Unknown workunit received. " << std::endl;
		}

	}

	save_state_auto();
	print_stats();
}



void
MPI_LoopHashRefine_Emperor::process_outbound_wus(){
}


// iterate through the structure store and add strucutres to the central library according to the algorithm specified or the default algorithm
bool
MPI_LoopHashRefine_Emperor::add_structures_to_library( SilentStructStore &new_structs, std::string add_algorithm ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool result = false;

	for( SilentStructStore::const_iterator it = new_structs.begin();
		 it != new_structs.end(); ++it )
	{
		runtime_assert( *it );
		core::io::silent::SilentStruct *pss = &(*(*it));

		// Filter for max_emperor_lib_round_

		if( max_emperor_lib_round_ > 0 ){
			core::Size structure_round = (core::Size) pss->get_energy("round");
			if( structure_round > max_emperor_lib_round_ ) continue;
		}

		// add the structure if it passes energy and rms filters evaluated further down there
		runtime_assert( pss );
		pss->add_energy( "emperor_count",  pss->get_energy( "emperor_count" ) );
		bool local_result = add_structure_to_library( *pss, add_algorithm );
		result |= local_result;
	}

	// always dump *everything* returned to emperor!
	dump_structures( new_structs, false );

	limit_library();
	print_library();
	return result;
}


} // namespace loophash
} // namespace protocols


