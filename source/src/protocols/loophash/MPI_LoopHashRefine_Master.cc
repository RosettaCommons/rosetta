// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Master.cc
/// @brief
/// @author Mike Tyka

#define TRDEBUG TR.Debug

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/loophash/MPI_LoopHashRefine.hh>
#include <protocols/loophash/MPI_LoopHashRefine_Master.hh>
#include <protocols/loophash/WorkUnit_LoopHash.hh>
#include <protocols/relax/WorkUnit_BatchRelax.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <core/pose/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <ObjexxFCL/format.hh>
/// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

#ifndef _WIN32 // REQUIRED FOR WINDOWS
#endif

#include <fstream>
#include <utility/string_util.hh>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loophash {

using namespace protocols::wum;

static thread_local basic::Tracer TR( "MPI.LHR.Master" );


void
MPI_LoopHashRefine_Master::set_defaults(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	max_loophash_per_structure_ = option[ OptionKeys::lh::max_loophash_per_structure ]();
	batch_relax_chunks_         = option[ OptionKeys::lh::mpi_batch_relax_chunks ]();
	batch_relax_absolute_max_   = option[ OptionKeys::lh::mpi_batch_relax_absolute_max ]();
	outbound_wu_buffer_size_    = option[ OptionKeys::lh::mpi_outbound_wu_buffer_size ]();
	loophash_split_size_        = option[ OptionKeys::lh::mpi_loophash_split_size     ]();
	library_expiry_time_        = option[ OptionKeys::lh::library_expiry_time ]();
	expire_after_rounds_        = option[ OptionKeys::lh::expire_after_rounds ]();
	mpi_master_save_score_only_ = option[ OptionKeys::lh::mpi_master_save_score_only ]();
}


void
MPI_LoopHashRefine_Master::init(){
	// Are we resuming an old job ?
	if( mpi_resume() != "" ){
		TR << "Resuming job from IDENT:  " <<  mpi_resume() << std::endl;
		load_state( mpi_resume() );
	} else {
		load_structures_from_cmdline_into_library( max_lib_size() * master_rank() );
	}
	
	// sample_weight cannot be initialized until after structures are imported, so we can check size
		load_sample_weight();
	TR << "STARTLIB: " << std::endl;
	print_library();
}


void
MPI_LoopHashRefine_Master::go()
{
	// initialize master (this is a virtual functino call and this function is overloaded by the children of this class)
	TR << "Init Master: " << mpi_rank() << std::endl;
	init();

	TR << "Master Node: Waiting for job requests..." << std::endl;
	while(true){
		// process any incoming messages such as incoming
		TRDEBUG << "Master: processing msgs.." << std::endl;
		process_incoming_msgs();

		TRDEBUG << "Master: process incoming" << std::endl;
		process_inbound_wus();

		TRDEBUG << "Master: process outbound" << std::endl;
		process_outbound_wus();

		// ok, we've done all our work, now wait until we hear from our slaves
		process_incoming_msgs( true );

		print_stats_auto();
	}
}


/// @brief figure out what to do with incoming WUs.
/// Some will be returning WUs that need to be resent others will be finished and will need
/// reintegration into the library
void
MPI_LoopHashRefine_Master::process_inbound_wus(){
	using namespace protocols::loops;

	check_library_expiry_dates();
	TRDEBUG << "Finished checking library dates"<<std::endl;
	if( inbound().size() > 0 ){
		TRDEBUG << "Processing inbound WUs on master.." << std::endl;
	}
	while( inbound().size() > 0 )
	{
		WorkUnitBaseOP  next_wu =  inbound().pop_next();
		runtime_assert( next_wu != 0 );

		// skip returning waiting WUs
		if ( next_wu->get_wu_type() == "waitwu" ) continue;

		// Upcast to a StructureModifier WU
		WorkUnit_SilentStructStoreOP structure_wu = utility::pointer::dynamic_pointer_cast< protocols::wum::WorkUnit_SilentStructStore > ( next_wu );

		// If upcast was unsuccessful - warn and ignore.
		if ( structure_wu.get() == NULL ){
			TR << "Cannot save structural data for WU: " << std::endl;
			next_wu->print( TR );
			continue;
		}

		// Otherwise extract structures and figure out what to do with them
		TRDEBUG << "Saving decoy store.. " << std::endl;
		SilentStructStore &decoys = structure_wu->decoys();

		if ( structure_wu->get_wu_type() == "loophasher" ){
			totaltime_loophash() += structure_wu->get_run_time();
			TR << "LoopHash return: " << decoys.size() << " structs in " << structure_wu->get_run_time() << "s " << " frm " << structure_wu->last_received_from() << std::endl;
			// Add the node that returned to the blacklist of all WUs with the same ssid and start_ir
			for( WorkUnitQueue::iterator iter = outbound().begin(); iter != outbound().end(); iter++ ) {
					if( (*iter)->get_wu_type() == "loophasher" ) {
			/*			// Upcast to a StructureModifier WU
						// Why use dynamic cast when we're sure of type? (copied from above)
						// So slow!
					TR << "im here1" << std::endl;
						WorkUnit_SilentStructStore* i = dynamic_cast<  WorkUnit_SilentStructStore * > ( (WorkUnitBase*) (&(*(*iter))) );
					TR << "im here2" << std::endl;
							if( ( core::Size )(*i->decoys().begin())->get_energy("ssid") == (core::Size)(*decoys.begin())->get_energy("ssid") ) {
//	i->extra_data_1() == structure_wu->extra_data_1() ) {
//	*/
								if( (*iter)->extra_data_1() == structure_wu->extra_data_1() && (*iter)->extra_data_3() == structure_wu->extra_data_3() ) {
									(*iter)->add_blacklist( structure_wu->last_received_from() );
									TRDEBUG << "Added node " << structure_wu->last_received_from() << " to blacklist of WU " << (*iter)->id() << std::endl;
							}
				  }
		  }
			n_loophash_++;
			//to_be_relaxed_.add( decoys );
			if( decoys.size() > 0 ){
				add_relax_batch( decoys );
				total_structures_ += decoys.size();
			}
		} else
		if ( structure_wu->get_wu_type() == "resultpack" ){
			decoys.all_sort_silent_scores();
			// dump structures
			TR << "Emperor sent: " << decoys.size() << " structs" << std::endl;
			print_library();
			add_structures_to_library( decoys, "add_n_limit" );
			print_library();
			// dont dump structures that came straight from emperor. 
			//dump_structures( decoys, mpi_master_save_score_only_ );
		} else
		if ( structure_wu->get_wu_type() == "batchrelax" ){
			decoys.all_sort_silent_scores();
			decoys.all_add_energy("state", 2 ); // mark structures are just having come thoruhg batchrelax
			totaltime_batchrelax_ += structure_wu->get_run_time();
			n_batchrelax_ ++;
			TR << "BatchRelax return: " << decoys.size() << " structs in " << structure_wu->get_run_time() << "s " << " frm " << structure_wu->last_received_from() << std::endl;
			add_structures_to_library( decoys );
			dump_structures( decoys, mpi_master_save_score_only_ );
		} else {
			TR.Error << "Unknown workunit received. " << std::endl;
		}


	}

	print_stats();
}


void
MPI_LoopHashRefine_Master::process_outbound_wus(){
	TRDEBUG << "Adding loophash WUs if necessary .. " << std::endl;
	if( outbound().size() < outbound_wu_buffer_size_ ){
		if ( library_central().size() == 0 ){
			TR.Error << "FATAL ERROR:  library_central_ is empty! " << std::endl;
      utility_exit_with_message( "FATAL ERROR:  library_central_ is empty! " );
		}
		// pick a random structure from the library

		core::Size finished_structures=0;
		for( SilentStructStore::iterator it = library_central().begin(); it !=  library_central().end(); it ++ ){
			if( max_loophash_per_structure_ > (*it)->get_energy("lhcount"))
			{
					TRDEBUG << "Adding: " << (*it) << "  " << (*it)->get_energy("lhcount") << std::endl;
				 (*it)->add_energy( "lhcount",  (*it)->get_energy("lhcount") + 1.0 );
				create_loophash_WUs( *it );
			}else{
					finished_structures += 1;
					TRDEBUG << "Already done: " << (*it) << "  " << (*it)->get_energy("lhcount") << std::endl;
			}
		}
		TR << "WARNING: " << finished_structures << "  " << library_central().size() << std::endl;
		if ( finished_structures >= library_central().size() ){
			TR << "WARNING: The starting structs exhausted!" << std::endl;
		}
	}

	save_state_auto();
}


void
MPI_LoopHashRefine_Master::create_loophash_WUs( const core::io::silent::SilentStructOP &start_struct ){

		runtime_assert( start_struct != 0 );
		core::pose::Pose start_pose;
		start_struct->fill_pose( start_pose );
  	core::util::switch_to_residue_type_set( start_pose, core::chemical::CENTROID);
		core::pose::set_ss_from_phipsi( start_pose );

		//refresh the sampling weight comment, as it may have changed
		// easier to do it here, copy_scores copies comments as well
		core::pose::delete_comment(start_pose,"sample_weight");
		core::pose::add_comment(start_pose,"sample_weight", sample_weight_str_);

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		
		core::io::silent::SilentStructOP ss = option[ OptionKeys::lh::bss]() ?
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary") :
				core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		ss->fill_struct( start_pose );
		ss->copy_scores( *start_struct );

		// first cound up "round" counter - just counts how many times each structure has been
		// thorugh the loop hasher
		core::Size round = (core::Size)ss->get_energy("round");
		round++;
		ss->add_energy("round", round );
		ss->add_energy("masterid", mpi_rank() );
		ss->add_energy("parent_score", ss->get_energy("score") );

		core::Size start_ir = 1;
		core::Size end_ir = 1;
		core::Size ssid = (core::Size)ss->get_energy("ssid");

		core::Size count_wus = 0;
		for( ;start_ir< start_pose.total_residue(); start_ir+=loophash_split_size_ )
		{
			end_ir =  std::min( start_ir + loophash_split_size_ - 1, start_pose.total_residue());
			if( end_ir < start_ir) end_ir = start_ir;
		  if( start_pose.total_residue() - end_ir < loophash_split_size_ ) end_ir = start_pose.total_residue();
			TRDEBUG << "Adding a new loophash WU: " << start_ir << " - " << end_ir << ", ssid = " << ssid << std::endl;
			
			count_wus++;
			WorkUnit_LoopHashOP new_wu( new WorkUnit_LoopHash( start_ir, end_ir, ssid ) );
			// this is unsatisfying.. why can't i use the template ? grrr C++ thou are limited.
			new_wu->set_wu_type("loophasher");
			new_wu->decoys().add( ss);
			new_wu->clear_serial_data();
			outbound().add( new_wu );
		  if( start_pose.total_residue() - end_ir < loophash_split_size_ ) start_ir = start_pose.total_residue();
		}
		TR << "Added " << count_wus << " loophash WUs to queue. ssid=" << ssid << std::endl;

}


void
MPI_LoopHashRefine_Master::add_relax_batch( SilentStructStore &start_decoys ){
	if( start_decoys.size() == 0 ) return;
	TR << "Adding relax WUs.." << start_decoys.size() << std::endl;

	core::Size count_adds = 0;
	core::Size count_adds_b4_limit = 0;
	core::Size count_wus = 0;

	core::Size chunks  = 1 + core::Size( floor( core::Real(start_decoys.size()) / core::Real( batch_relax_chunks_ ) ) );
	core::Size batchrelax_batchsize_ = (start_decoys.size() / chunks) + 1;
	core::Size dcount=0;
	while( dcount < start_decoys.size() ){
		protocols::relax::WorkUnit_BatchRelaxOP new_wu( new protocols::relax::WorkUnit_BatchRelax_and_PostRescore() );
		new_wu->set_wu_type("batchrelax");
		core::Size lcount=0;

		for(lcount=0; lcount < batchrelax_batchsize_; lcount++ ){
			if ( dcount < start_decoys.size() ){
				core::io::silent::SilentStructOP new_relax_structure =  start_decoys.get_struct( dcount );
				TRDEBUG << "AddRelaxStructure: " << format_silent_struct(new_relax_structure)  << std::endl;
				new_wu->decoys().add( new_relax_structure );
			}
			dcount++;
		}

		// Mix up the order
    //std::random__shuffle( new_wu->decoys().begin(), new_wu->decoys().end());
    numeric::random::random_permutation(new_wu->decoys().begin(), new_wu->decoys().end(), numeric::random::rg());

		// make sure the chunk size doesnt exceed batch_relax_absolute_max_
		core::Size chunk_size = new_wu->decoys().size();
		new_wu->decoys().limit( batch_relax_absolute_max_ );

		total_structures_relax_ += new_wu->decoys().size();
		new_wu->clear_serial_data();

		count_adds += new_wu->decoys().size();
		count_adds_b4_limit += chunk_size;
		count_wus ++;
		// Relax work units have a lot of structures and fill up the queue and lead to memory crashes. Thus they get prioritized and added at the beginning of the queue!!
		outbound().push_front( new_wu );
	}

	TR << "Adding " << count_adds << "/" << count_adds_b4_limit << " structs for batchrlx. " << count_wus << " WUs" << std::endl;

}


// this goes through the library and identifies structures that have not managed to get replaced
// for some cutoff amount of time. It will send back this structure and request a new structure with the same ssid from
// the emperor.
void
MPI_LoopHashRefine_Master::check_library_expiry_dates(){
	core::Size current_time = time(NULL);

	SilentStructStore::iterator jt_last = library_central().begin();

	for( SilentStructStore::iterator jt =  library_central().begin();
									jt != library_central().end(); jt ++ )
	{
		TR.Debug << "Checking structure.." << std::endl;
		core::Size struct_time = (core::Size)(*jt)->get_energy("ltime");
		core::Size ssid        = (core::Size)(*jt)->get_energy("ssid");
		core::Size round       = (core::Size)(*jt)->get_energy("round");
		
		bool expired = false;
		// is the structure expired due to time limit ?		
		if( (int(current_time) - int(struct_time)) > (int)library_expiry_time_ ){
			expired = true;
			TR << "Structure: " << ssid << " is expired: " << int(current_time) - int(struct_time) << " > " << (int)library_expiry_time_ <<  std::endl;
		}

		// is the structure expired because it has done too many rounds ?
		if( (expire_after_rounds_ > 0) && ( round >= expire_after_rounds_ ) ){
			expired = true; 
			TR << "Structure: " << ssid << " Round:  is expired: " << round << " >= " << expire_after_rounds_ << std::endl; 
		}

		if( ! expired ){
			jt_last = jt;
			continue;
		}

		// ok, so the structure is expired. send it to the emperor and kill it. wait for a new structure to arrive

		(*jt)->add_energy("expire", (core::Size)(*jt)->get_energy("expire") + 1);
		

		// send the expired structure to the emperor (who will in due time send back a new one)
		WorkUnit_SilentStructStoreOP getnewstruct( new WorkUnit_SilentStructStore( ) );
		getnewstruct->set_wu_type( "getnewstruct" );
		getnewstruct->decoys().add( (*jt) );
		send_MPI_workunit( getnewstruct, 0 ); // The 0 is the MPI_RANK of the master - constant would be better here!
		
		// clear the queue of loophash WUs from previous struct to avoid false blacklisting
		// assume that false blacklisting from currently processing loophash WU is unlikely 
		
		core::Size erase_count = 0;
		for( WorkUnitQueue::iterator iter = outbound().begin(); iter != outbound().end();) {
				if( (*iter)->get_wu_type() == "loophasher" && ssid == (*iter)->extra_data_3() ) {
					TRDEBUG<<"erasing wu" <<std::endl;
					iter->reset(); // to NULL
					TRDEBUG<<"erasing wu from list" <<std::endl;
					iter = outbound().erase( iter );
					TRDEBUG<<"erasing done" <<std::endl;
					erase_count ++;
				} else {
						++iter;
				}
		}
		TR << "Erased " << erase_count << " deprecated WUs from outbound queue" << std::endl;
	
		// now delete this expired structure - it is now at the emperor's mercy
		library_central().erase(jt);

		TR << "Reported expired structure to emperor: - waiting for new structure" << std::endl;
		receive_MPI_workunit( 0 ); //receive the reply from master and add it to the normal inbound queue. the 0 here is the emperor's MPIRANK. Better replace with a function or constant
		TR << "Done. Restarting reporting.." << std::endl;	
		break; // only one at a time.
		// reset the iterator to the beginning - we must do that because we could have added the new structure whereever - beginning is the only save iterator 
		jt=library_central().begin();
			
		TRDEBUG << "Library state: " << std::endl;	
		print_library();
	}
	TRDEBUG << "end of check_library_expiry_dates" << std::endl;
}

/// This is a virtual over load of the base class MPI_LoopHashRefine:: add_structure_to_library with an extra behavioural step
/// that reports any successful library add-ons to the emperor. This behaviour is master specific and thus should not be in the base class.

bool
MPI_LoopHashRefine_Master::add_structure_to_library( core::io::silent::SilentStruct &pss, std::string add_algorithm ){
	bool result = MPI_LoopHashRefine::add_structure_to_library( pss, add_algorithm );
	TR << "MPI_LoopHashRefine_Master::add_structure_to_library: " << std::endl;
	if(result) report_structure_to_emperor( pss );
	return result;
}

void
MPI_LoopHashRefine_Master::report_structure_to_emperor(  core::io::silent::SilentStructOP &ss ) {
	WorkUnit_SilentStructStoreOP resultpack( new WorkUnit_SilentStructStore( ) );
	resultpack->set_wu_type( "resultpack" );
	resultpack->decoys().add( ss );
	send_MPI_workunit( resultpack, my_emperor() );
	TR << "Reported structure to emperor: " << format_silent_struct( ss ) << std::endl;
}

void
MPI_LoopHashRefine_Master::report_structure_to_emperor(  core::io::silent::SilentStruct &pss ) {
	WorkUnit_SilentStructStoreOP resultpack( new WorkUnit_SilentStructStore( ) );
	resultpack->set_wu_type( "resultpack" );
	resultpack->decoys().add( pss );
	send_MPI_workunit( resultpack, my_emperor() );
	TR << "Reported structure to emperor: " << format_silent_struct(pss) << std::endl;
}


void
MPI_LoopHashRefine_Master::load_sample_weight() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// This just loads sample weights from a file
		// I assume that the optionkeys sanitizes input
		
		// using ifstream instead of utility::io::izstream because izstream doesn't return success bool 
		if( option[ OptionKeys::lh::sample_weight_file ].active() ) {
				std::string pathtofile = option[ OptionKeys::lh::sample_weight_file ]();
				std::ifstream file( pathtofile.c_str() ); 
				if (!file) utility_exit_with_message( "Failed to open sample_weight file.  Check path." );
				std::string line, tmp;
				while(getline( file, line ) ) {

					boost::trim(line);
					std::vector < std::string > r;
					boost::split(r, line, boost::is_any_of("\t "));
					
				
					core::Real i=0.0;
					// Check for correct format
					try {
						i = boost::lexical_cast<core::Real> (r[1] );
					} catch( boost::bad_lexical_cast &) {
							utility_exit_with_message( "Sample weight second column can't be casted to an int.");
					}

					if (i < 0) {
							utility_exit_with_message( "Sample weight second column is not an float larger than 0." );
						} else {
							tmp += r[1] + " ";
						}
				}
				// check for correct length
				boost::trim(tmp);
				std::list < std::string > t;
				t = utility::split_to_list(tmp);
				if( t.size() != (*(library_central().begin()))->nres() )
					utility_exit_with_message( "Sample weight file either improperly formatted or does not have same number of residues as structure." );
				TR << "Sample weight file successfully loaded" << std::endl;
				sample_weight_str_ = tmp;
		} else {
				TR << "Using default sample weight of 50 for every residue" << std::endl;
				std::string t = "50";
				for( Size i = 0; i < (*(library_central().begin()))->nres() - 1; i++ ) {
						t += " 50";
				}
				sample_weight_str_ = t;
		}
}


} // namespace loophash
} // namespace protocols


