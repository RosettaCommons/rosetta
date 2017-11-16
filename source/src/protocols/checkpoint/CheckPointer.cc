// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Oliver Lange
/// @author Mike Tyka

// Project headers
// This has to come before boinc_util.hh or we get this error on VC++
// '_read' : is not a member of 'std::basic_istream<_Elem,_Traits>'
#include <utility>
#include <utility/io/izstream.hh>

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#ifdef BOINC
#include <protocols/boinc/boinc.hh>
#ifndef ANDROID // do not use watchdog thread for ANDROID
#include <protocols/boinc/watchdog.hh>
#endif
#endif

// Unit Headers
#include <protocols/checkpoint/CheckPointer.hh>

// Package Headers
#include <protocols/checkpoint/Checkpoint.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_writer.hh>

#if defined(WIN32) || defined(BOINC)
#include <core/io/silent/SilentStructFactory.hh>
#endif

#include <core/io/silent/SilentStruct.hh>
#include <numeric/random/random.hh>

// ObjexxFCL Headers
// Utility headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>

#ifdef BOINC
#ifndef _WIN32
#include "pthread.h"
#endif
#endif


// C++ headers

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>



using basic::Error;
using basic::Warning;

namespace protocols {
namespace checkpoint {

static basic::Tracer TR( "protocols::checkpoint" );

using namespace core;


void FileBuffer::dump(){

	//std::cerr << "Dumped checkpoint. " << std::endl;
	if ( !gzipped_ ) {
		utility::io::ozstream output( filename_ );
		if ( !output ) {
			std::cout << "Cannot write checkpoint! open failed for file: " << filename_ << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		output << contents_;
		output.close();
	} else {
		utility::io::ozstream output( filename_ );
		output << contents_;
		output.close();
	}

}


bool pose_to_binary_silent_file( std::ostream &output, const std::string &tag, const pose::Pose &pose ){
	using namespace io::silent;
	SilentFileOptions opts;
	SilentFileData outsfd( opts );
	SilentStructOP pss( new BinarySilentStruct( opts ) );
	pss->fill_struct( pose, tag );

	pss->print_header( output );
	pss->print_scores( output );
	pss->print_conformation( output );

	return true;
}

bool pose_from_binary_silent_file( const std::string &filename, const std::string &tag, pose::Pose &pose, bool fullatom ){
	using namespace core::chemical;
	ResidueTypeSetCOP residue_set;
	pose::Pose tmppose;
	if ( fullatom ) residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	else           residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );
	core::io::silent::SilentFileOptions opts;
	core::io::silent::SilentFileData sfd("", false, false, "binary", opts );
	if ( !sfd.read_file( filename ) ) return false;
	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		if ( iter->decoy_tag() != tag ) continue;
		iter->fill_pose( tmppose, *residue_set );
		break;
	}
	tmppose.transfer_constraint_set( pose );
	pose.clear();
	pose = tmppose;
	return true;
}

CheckPointer::CheckPointer( std::string  type ):
	type_ (std::move( type )),
	disabled_( false ),
	count_checkpoint_recoveries_( 0 )
{
	using namespace basic::options;
	delete_checkpoints_ = option[ OptionKeys::run::delete_checkpoints ]();
	if ( option[ OptionKeys::run::suppress_checkpoints ]() ) disabled_ = true;
}


void CheckPointer::debug( const std::string &tag, const std::string &label, core::Real data1, core::Real data2, core::Real data3 ) const {
	using namespace basic::options;
	if ( disabled_ ) return;
	if ( option[ OptionKeys::run::checkpoint ]() ||
			option[ OptionKeys::run::checkpoint_interval ].user() ) {

		TR.Info << "CHECKDEBUG: " << tag << "_" << label << "  " << data1 << "  " << data2 << "  " << data3 << std::endl;
	}
}

void CheckPointer::flush_checkpoints()
{
	for ( auto & i : file_buffer ) {
		i.dump();
	}
	file_buffer.clear();
}

core::Size CheckPointer::file_buffer_size()
{
	core::Size total_size = 0;
	for ( auto & i : file_buffer ) {
		total_size += i.size();
	}
	return total_size;
}

void CheckPointer::checkpoint(
	pose::Pose &pose,
	moves::MonteCarlo* mc,
	std::string const& current_tag,
	std::string const& id,
	bool /*foldtree*/
)
{


	///////////// Is checkpointing on at all ?

	if ( disabled_ ) return;
	if ( !Timer::is_on() ) return;

	if ( pose.size() <= 0 ) return;

#ifdef BOINC
#ifndef ANDROID

	/////////////
	// set global bailout structure - this is used by the watchdog to abort when the watchdog time
	// has been  global_bailout = new pose::Pose( pose );


	boinc_begin_critical_section();

	using namespace core::io::silent;

	std::stringstream bail_structure;
	std::stringstream bail_structure_header;
	SilentFileOptions opts;
	SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct_out( opts );
	ss->fill_struct( pose,  "W_00000001" );
	SilentFileDataOP sfd( new SilentFileData( opts ) );

	// write the structure to the bailout buffer
	sfd->write_silent_struct( *ss, bail_structure );

	// write the header to the bailout header bufffer (might be needed later if the watchdog
	// is triggered before the first decoy is written.

	ss->print_header( bail_structure_header );

	// make sure the watchdog isn't about to use that variable!
#ifdef _WIN32
#else
	pthread_mutex_lock(&protocols::boinc::watchdog::bailout_mutex);
#endif
	protocols::boinc::watchdog::bailout_silent_structure = bail_structure.str();
	protocols::boinc::watchdog::bailout_silent_structure_header = bail_structure_header.str();
#ifdef _WIN32
#else
	pthread_mutex_unlock(&protocols::boinc::watchdog::bailout_mutex);
#endif

	boinc_end_critical_section();
#endif // ANDROID
#endif // BOINC

	////// Create Checkpoint in buffer.


	/// Create the unique checkpoint ID
	std::string mcstr( "_MC" );
	if ( mc == nullptr ) mcstr = "";
	std::string checkpoint_id( "chk_" + current_tag + "_" + type() + "_" + mcstr + "_" + id );
	if ( pose.is_fullatom() ) checkpoint_id += "_fa";

	/// Make a note of it.
	checkpoint_ids_.push_back( checkpoint_id );

	/// Save Random State

	{
		FileBuffer new_file( checkpoint_id + ".rng.state.gz", true /*gzipped*/ );
		std::stringstream ss_stream;
		numeric::random::rg().saveState( ss_stream );
		new_file.set_contents( ss_stream.str() );
		file_buffer.push_back( new_file );
	}
	//  utility::io::ozstream ozs(checkpoint_id + ".rng.state.gz");
	//  numeric::random::RandomGenerator::saveAllStates(ozs);
	//  ozs.close();


	//std::string notag="";
	if ( mc != nullptr ) {
		{
			FileBuffer new_file(  checkpoint_id + "mc_last.out"  );
			std::stringstream ss_stream;
			pose_to_binary_silent_file( ss_stream , checkpoint_id, mc->last_accepted_pose() );
			new_file.set_contents( ss_stream.str() );
			file_buffer.push_back( new_file );
		}

		{
			FileBuffer new_file(  checkpoint_id + "mc_low.out"  );
			std::stringstream ss_stream;
			pose_to_binary_silent_file( ss_stream , checkpoint_id, mc->lowest_score_pose() );
			new_file.set_contents( ss_stream.str() );
			file_buffer.push_back( new_file );
		}
		//  pose_to_binary_silent_file( checkpoint_id + "mc_last.out", checkpoint_id, mc->last_accepted_pose() );
		//  pose_to_binary_silent_file( checkpoint_id + "mc_low.out",  checkpoint_id, mc->lowest_score_pose() );
	}

	{
		FileBuffer new_file(  checkpoint_id + ".out"  );
		std::stringstream ss_stream;
		pose_to_binary_silent_file( ss_stream , checkpoint_id, pose );
		new_file.set_contents( ss_stream.str() );
		file_buffer.push_back( new_file );
	}


	TR << "CHECK: " << " Created virtual checkpoint: " << checkpoint_id << std::endl;


	/////////////////////////////////////////////////////////////////////////////////////////
	//// Is it time to flush the checkpoint ?


	if ( protocols::checkpoint::Timer::time_to_checkpoint() ) {
		// ------ BOINC critical ------------------------------------
#ifdef BOINC
			boinc_begin_critical_section();
#endif

		flush_checkpoints();

		// ------ endof BOINC critical (Timer::reset calls boinc_end_critical_section() )------------------------------------
		protocols::checkpoint::Timer::reset(); // required for checkpoint
	} else {

		// WHatever the user says - never store more than 10 checkpoints - it takes too much memory. Allow 15MB
		// but no more. if we exceed that then dump checkpoint contents whatever.

		if ( file_buffer_size() > 15000000 ) {
			flush_checkpoints();
		}

	};

	//std::cerr << "CHECK: CheckPointerSize: " << file_buffer_size() << "   N: " << file_buffer.size() << std::endl;
}

bool CheckPointer::recover_checkpoint(
	pose::Pose & pose, moves::MonteCarlo* mc,
	std::string const & current_tag,
	std::string const & id,
	bool fullatom,
	bool foldtree
)
{
	using namespace basic::options;
	using namespace core::chemical;

	bool debug = option[ OptionKeys::run::debug ]();

	if ( disabled_ ) return false;
	if ( !protocols::checkpoint::Timer::is_on() ) return false; // required for checkpoint

	TR.Error  << "recovering checkpoint of tag " << current_tag
		<< " with id " << id
		<< std::endl;

	std::string mcstr( "_MC" );
	if ( mc == nullptr ) mcstr = "";
	// must be same as in checkpoint function
	std::string checkpoint_id( "chk_" + current_tag + "_" + type() + "_" + mcstr + "_" + id );
	if ( fullatom ) checkpoint_id += "_fa";


	ResidueTypeSetCOP residue_set;


	// Make sure the check point exists - presence of the silent file indicates so.

	TR << "CHECKING FOR CHECKPOINT" <<  " Id: " << checkpoint_id << " : ";
	if ( !utility::file::file_exists( checkpoint_id + ".out" ) ) {
		TR << " NOT PRESENT" << std::endl;
		return false;
	}

	// Make an appropriate residue set.
	if ( fullatom ) residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	else           residue_set = ChemicalManager::get_instance()->residue_type_set( CENTROID );

	// Restore random number state!
	if ( utility::file::file_exists( checkpoint_id + ".rng.state.gz" ) ) {
		TR.Info << "Read random number state. " << std::endl;
		utility::io::izstream izs(checkpoint_id + ".rng.state.gz");
		numeric::random::rg().restoreState(izs);
		izs.close();
	} else {
		utility_exit_with_message("Random generator state not found for checkpoint " + checkpoint_id + " even though the silent file is there. " );
	}

	// get the structure from the binary file
	pose_from_binary_silent_file( checkpoint_id + ".out", checkpoint_id, pose, fullatom );

	if ( pose.is_fullatom() != fullatom ) {
		utility_exit_with_message("Fullatom mismatch in checkpointer.");
	}

	if ( debug ) {
		core::io::StructFileRepOptionsOP sfr_opts( new core::io::StructFileRepOptions() );
		sfr_opts->set_fold_tree_io( foldtree );
		core::io::pdb::dump_pdb( pose, checkpoint_id + ".debug.pdb", sfr_opts );

	}

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif
	if ( mc != nullptr && utility::file::file_exists( checkpoint_id + ".mc_last.pdb" ) ) {
		pose::Pose recovered_mc_last =  mc->last_accepted_pose();
		pose_from_binary_silent_file( checkpoint_id + ".out", checkpoint_id, recovered_mc_last, fullatom );
		mc->set_last_accepted_pose( recovered_mc_last );
		if ( debug ) {
			core::io::StructFileRepOptionsOP sfr_opts( new core::io::StructFileRepOptions() );
			sfr_opts->set_fold_tree_io( foldtree );
			core::io::pdb::dump_pdb( mc->last_accepted_pose(), checkpoint_id + ".mc_last.debug.pdb"  );
		}
	}
	if ( mc != nullptr && utility::file::file_exists( checkpoint_id + ".mc_low.pdb" ) ) {
		pose::Pose recovered_mc_low =  mc->lowest_score_pose();
		pose_from_binary_silent_file( checkpoint_id + ".out", checkpoint_id, recovered_mc_low, fullatom );
		mc->set_lowest_score_pose( recovered_mc_low );
		if ( debug ) {
			core::io::StructFileRepOptionsOP sfr_opts( new core::io::StructFileRepOptions() );
			sfr_opts->set_fold_tree_io( foldtree );
			core::io::pdb::dump_pdb( mc->lowest_score_pose(), checkpoint_id + ".mc_low.debug.pdb" );
		}
	}
	checkpoint_ids_.push_back( checkpoint_id );
	TR << " SUCCESS" << std::endl;

	count_checkpoint_recoveries_ += 1;

	return true;
}

void CheckPointer::clear_checkpoints() {

	TR.Info << "Deleting checkpoints of " << type() << std::endl;
	if ( disabled_ ) return;
	using namespace basic::options;
	if ( delete_checkpoints_ ) {
		for ( auto & checkpoint_id : checkpoint_ids_ ) {
			//std::cerr << "deleting checkpoint files with id: " << checkpoint_ids_[i] << std::endl;
			utility::file::file_delete( checkpoint_id + ".mc_last.out" );
			utility::file::file_delete( checkpoint_id + ".mc_low.out" );
			utility::file::file_delete( checkpoint_id + ".out" );
			utility::file::file_delete( checkpoint_id + ".rng.state.gz" );
		}
	} else {
		TR.Debug << "Checkpoint deletion disabled!" << std::endl;
	}

	file_buffer.clear(); // also make sure no more structures get written out accidentally.
	checkpoint_ids_.clear();
} // ClassicAbinitio::clear_checkpoints()


}
}
