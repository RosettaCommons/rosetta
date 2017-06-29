// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ResidueTypeSetTests.cxxtest.hh
/// @brief unit test for ResidueTypeFinder
/// @author Rocco Moretti (rmorettiase@gmail.com),

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueTypeSetCache.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>
#include <ostream>
#include <vector>

#ifdef MULTI_THREADED

// Utility headers
#include <utility/thread/ReadWriteMutex.hh>

// C++11 headers
#include <thread>
#include <chrono>

#endif

static basic::Tracer TR("core.chemical.ResidueTypeCacheConcurrencyTests.cxxtest");


using namespace core::chemical;

class DummyGRTS : public GlobalResidueTypeSet
{
public:
	DummyGRTS(
		std::string const & name,
		std::string const & directory
	) : GlobalResidueTypeSet( name, directory ) {}

	virtual ~DummyGRTS() {}

#ifdef MULTI_THREADED
	void pick_up_read_lock_on_cache() {  cache_object()->read_write_mutex().obtain_read_lock(); }
	void release_read_lock_on_cache() {  cache_object()->read_write_mutex().release_read_lock(); }
	void pick_up_write_lock_on_cache() { cache_object()->read_write_mutex().obtain_write_lock(); }
	void release_write_lock_on_cache() { cache_object()->read_write_mutex().release_write_lock(); }
#endif
};

typedef utility::pointer::shared_ptr< DummyGRTS > DummyGRTSOP;

class ResidueTypeSetCacheConcurrencyTests : public CxxTest::TestSuite {

	DummyGRTSOP rts_;
	bool finished_;
	bool started_;

public:

	void setUp() {
		finished_ = false;
		started_ = false;
		core_init();
#ifdef MULTI_THREADED
		// only bother loading the RTS if we're actually going to use it
		rts_ = DummyGRTSOP( new DummyGRTS( "fa_standard", basic::database::full_name( "chemical/residue_type_sets/fa_standard/" ) ));

#endif

	}

	void tearDown() {
		// deallocate the ResidueTypeSet
		rts_ = 0;
	}

	/// @brief Load a residue type from the RTS -- if the ALA:CtermProteinFull RT has not
	/// been created yet for the RTS, then this will obtain a write lock
	void
	ask_for_ala_cterm_full()
	{
		started_ = true;

		ResidueTypeCOP ala = rts_->name_mapOP( "ALA:CtermProteinFull" );

		// this assertion is not to test concurrency -- simply that we've set up
		// the RTS correctly; we should find this patched residue type
		TS_ASSERT( ala );

		finished_ = true;
	}

	/// @brief Load a residue type from the RTS -- if the PHE:CtermProteinFull RT has not
	/// been created yet for the RTS, then this will obtain a write lock
	void
	ask_for_phe_cterm_full()
	{
		started_ = true;

		ResidueTypeCOP phe = rts_->name_mapOP( "PHE:CtermProteinFull" );

		// this assertion is not to test concurrency -- simply that we've set up
		// the RTS correctly; we should find this patched residue type
		TS_ASSERT( phe );

		finished_ = true;
	}

	/// @brief Load a residue type from the RTS.
	void
	find_restypes_w_restype_finder_1()
	{
		started_ = true;

		//std::cout << "asking for ala:cterm" << std::endl;
		ResidueTypeCOP ala = rts_->name_mapOP( "ALA:CtermProteinFull" );
		//std::cout << "asking for all types with variants aa_phe" << std::endl;
		try {
			// internally, the RTS uses a ResidueTypeFinder in this method.
			rts_->get_all_types_with_variants_aa( aa_phe, ala->variant_types() );
		} catch ( utility::excn::EXCN_Base & e ) {
			std::cerr << "Caught exception " << e.msg() << std::endl;
		}
		finished_ = true;
	}

	void
	can_ser_turn_into_sep()
	{
		started_ = true;
		rts_->generates_patched_residue_type_with_name3( "SER", "SEP" );
		finished_ = true;
	}

	void
	can_cys_turn_into_scy()
	{
		started_ = true;
		rts_->generates_patched_residue_type_with_interchangeability_group( "CYS", "SCY" );
		finished_ = true;
	}

	void wait_for_success_and_kill_on_deadlock( std::string const & funcname )
	{
#ifdef MULTI_THREADED
		// After five seconds of the finished_ variable not being set to true,
		// declare deadlock, and exit.

		for ( int ii = 1; ii <= 10; ++ii ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(500) );
			if ( finished_ ) break;
		}
#endif

		TS_ASSERT( finished_ );
		if ( ! finished_ ) {
			std::cerr << "ERROR! ResidueTypeSetCacheConcurrencyTests::" << funcname << " deadlocked!" << std::endl;
		}
	}


	void test_load_new_restype_single_thread() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type
		// and wait a second; if that thread hasn't returned, it's
		// deadlocked and the test has failed.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::ask_for_ala_cterm_full, this );
		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_single_thread" );
		if ( finished_ ) running_thread.join();
#endif

	}

	void test_restype_finder_loading_new_restypes_in_single_thread() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		//std::cout << "Beginning test_restype_finder_loading_new_restypes_in_single_thread()" << std::endl;


		// Using the ResidueTypeFinder, look for ResidueTypes that are not already
		// in the RTS.
		// launch a thread that tries to get the new residue type
		// and wait a second; if that thread hasn't returned, it's
		// deadlocked and the test has failed.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		wait_for_success_and_kill_on_deadlock( "test_restype_finder_loading_new_restypes_in_single_thread" );
		if ( finished_ ) running_thread.join();

#endif

	}

	void test_load_new_restype_after_read_lock_released_two_threads_name_mapOP() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type.
		// That thread should begin running, but should block once it
		// tries to load the ala-cterm RT, and will have to wait until
		// the read lock is released.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::ask_for_ala_cterm_full, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// The other thread should not have been able to create the ala-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_read_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_after_read_lock_released_two_threads_name_mapOP" );
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_load_new_restype_after_write_lock_released_two_threads_name_mapOP() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		rts_->pick_up_write_lock_on_cache();

		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type.
		// That thread should begin running, but should block once it
		// tries to load the ala-cterm RT, and will have to wait until
		// the read lock is released.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::ask_for_ala_cterm_full, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// The other thread should not have been able to create the ala-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_write_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_after_write_lock_released_two_threads_name_mapOP" );
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_ser_to_sep_query_after_read_lock_released_two_threads() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		rts_->pick_up_read_lock_on_cache();

		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::can_ser_turn_into_sep, this );

		// wait for a bit for the thread to get started
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// The other thread should not have been able to finish until the read lock is released
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_read_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_ser_to_sep_query_after_read_lock_released_two_threads" );
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_cys_to_scy_query_after_read_lock_released_two_threads() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		rts_->pick_up_read_lock_on_cache();

		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::can_cys_turn_into_scy, this );

		// wait for a bit for the thread to get started
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// sleep a little while longer to give the other thread a chance to finish
		// if it were able (it shouldn't be able)
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );

		// The other thread should not have been able to finish until the read lock is released
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_read_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_cys_to_scy_query_after_read_lock_released_two_threads" );
		if ( finished_ ) running_thread.join();

#endif
	}


	void test_two_treads_read_from_RTS_name_mapOP() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create the ala-cterm RT
		ask_for_ala_cterm_full();
		started_ = finished_ = false; // reset progress booleans

		// Now create a read lock on the RTS
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that should now be in the RTS
		// launch a thread that tries to get the new residue type
		// That other thread may very well complete immediately
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::ask_for_ala_cterm_full, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}
		// ok, the other thread has gotten started; give it enough time to complete.
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );

		wait_for_success_and_kill_on_deadlock( "test_two_treads_read_from_RTS_name_mapOP" );
		if ( finished_ ) running_thread.join();
#endif
	}

	void test_two_treads_read_from_RTS_gen_patch_ser_to_sep() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		can_ser_turn_into_sep();
		started_ = finished_ = false; // reset progress booleans

		// Now create a read lock on the RTS
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that should now be in the RTS
		// launch a thread that tries to get the new residue type
		// That other thread may very well complete immediately
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::can_ser_turn_into_sep, this );

		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}
		// ok, the other thread has gotten started; give it enough time to complete.
		wait_for_success_and_kill_on_deadlock( "test_two_treads_read_from_RTS_gen_patch_ser_to_sep" );

		rts_->release_read_lock_on_cache();
		if ( finished_ ) running_thread.join();
#endif
	}

	void test_two_treads_read_from_RTS_gen_patch_cys_to_scy() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		can_cys_turn_into_scy();
		started_ = finished_ = false; // reset progress booleans

		// Now create a read lock on the RTS
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that should now be in the RTS
		// launch a thread that tries to get the new residue type
		// That other thread may very well complete immediately
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::can_cys_turn_into_scy, this );

		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}
		// ok, the other thread has gotten started; give it enough time to complete.
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );

		// The other thread should have been able to read from the RTS even while this thread
		// holds a read lock.
		wait_for_success_and_kill_on_deadlock( "test_two_treads_read_from_RTS_gen_patch_cys_to_scy" );

		rts_->release_read_lock_on_cache();
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_load_new_restype_after_read_lock_released_two_threads_restypefinder() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create the ala-cterm RT so that call does not block when we ask for the PHE:CtermFull RT
		ask_for_ala_cterm_full();
		started_ = finished_ = false; // reset progress booleans

		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type.
		// That thread should begin running, but should block once it
		// tries to load the phe-cterm RT, and will have to wait until
		// the read lock is released.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		std::this_thread::sleep_for( std::chrono::milliseconds(100) );

		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_read_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_after_read_lock_released_two_threads_restypefinder" );
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_load_new_restype_after_read_lock_released_five_threads_restypefinder() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create the ala-cterm RT so that call does not block when we ask for the PHE:CtermFull RT
		ask_for_ala_cterm_full();
		started_ = finished_ = false; // reset progress booleans

		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type.
		// That thread should begin running, but should block once it
		// tries to load the phe-cterm RT, and will have to wait until
		// the read lock is released.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		std::this_thread::sleep_for( std::chrono::milliseconds(100) );

		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Now lets release one read lock, there still should be three read locks remaining
		rts_->release_read_lock_on_cache();
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );
		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Now lets release another read lock, there still should be two read locks remaining
		rts_->release_read_lock_on_cache();
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );
		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Now lets release another read lock, there still should be one read lock remaining
		rts_->release_read_lock_on_cache();
		std::this_thread::sleep_for( std::chrono::milliseconds(100) );
		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Ok, release the final read lock to let the other thread get going.
		rts_->release_read_lock_on_cache();

		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_after_read_lock_released_five_threads_restypefinder" );
		if ( finished_ ) running_thread.join();

#endif
	}

	void test_load_new_restype_after_write_lock_released_two_threads_restypefinder() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create the ala-cterm RT so that call does not block when we ask for the PHE:CtermFull RT
		ask_for_ala_cterm_full();
		started_ = finished_ = false; // reset progress booleans

		rts_->pick_up_write_lock_on_cache();

		// ask for a ResidueType that's not already in the RTS
		// launch a thread that tries to get the new residue type.
		// That thread should begin running, but should block once it
		// tries to load the phe-cterm RT, and will have to wait until
		// the read lock is released.
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		// wait for a bit -- in the mean time, the ask_for_ala thread should have woken up,
		// then tried to obtain a write lock, and been unable to get one, and should be waiting
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		// The other thread should not have been able to create the phe-cterm RT as there is a read lock
		// held on the RTSC.
		TS_ASSERT( ! finished_ );

		// Ok, release the read lock to let the other thread get going.
		rts_->release_write_lock_on_cache();

		// wait for longer than you should have to so the other thread has an opportunity to run.
		std::this_thread::sleep_for( std::chrono::seconds(1) );

		wait_for_success_and_kill_on_deadlock( "test_load_new_restype_after_write_lock_released_two_threads_restypefinder" );
		if ( finished_ ) running_thread.join();

#endif
	}


	void test_two_treads_read_from_RTS_restypefinder() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create both ala-cterm- and then a whole bunch of c-term RTs for PHE
		// so the RTF won't cause the RTS to need any write locks
		find_restypes_w_restype_finder_1();
		started_ = finished_ = false; // reset progress booleans

		// Now create a read lock on the RTS
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that should now be in the RTS
		// launch a thread that tries to get the new residue type
		// That other thread may very well complete immediately
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		// wait for a bit for the other thread to have woken up,
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// Ok, the other thread should quickly complete even while the read lock remains in place.
		wait_for_success_and_kill_on_deadlock( "test_two_treads_read_from_RTS_restypefinder" );

		rts_->release_read_lock_on_cache();
		if ( finished_ ) running_thread.join();
#endif
	}

	void test_five_treads_read_from_RTS_restypefinder() {
		TS_ASSERT( true ); // for non-c++11-thread tests
#ifdef MULTI_THREADED
		// Create both ala-cterm- and then a whole bunch of c-term RTs for PHE
		// so the RTF won't cause the RTS to need any write locks
		find_restypes_w_restype_finder_1();
		started_ = finished_ = false; // reset progress booleans

		// Now create a bunch of read locks on the RTS
		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();
		rts_->pick_up_read_lock_on_cache();

		// ask for a ResidueType that should now be in the RTS
		// launch a thread that tries to get the new residue type
		// That other thread may very well complete immediately
		std::thread running_thread( &ResidueTypeSetCacheConcurrencyTests::find_restypes_w_restype_finder_1, this );

		// wait for a bit for the other thread to have woken up,
		while ( ! started_ ) {
			std::this_thread::sleep_for( std::chrono::milliseconds(1) );
		}

		// Ok, the other thread should quickly complete even while the read locks remains in place.
		wait_for_success_and_kill_on_deadlock( "test_five_treads_read_from_RTS_restypefinder" );

		rts_->release_read_lock_on_cache();
		rts_->release_read_lock_on_cache();
		rts_->release_read_lock_on_cache();
		rts_->release_read_lock_on_cache();
		if ( finished_ ) running_thread.join();
#endif
	}

};
