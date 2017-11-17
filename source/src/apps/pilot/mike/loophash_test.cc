// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief

// libRosetta headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


// C++ headers
#include <iostream>
// option key includes

#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#ifdef USEMPI
#include <mpi.h>
#endif

static basic::Tracer TR( "main" );


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


protocols::loophash::BackboneSegment get_backbone_segment(  protocols::loophash::LoopHashLibraryOP loop_hash_library,  numeric::geometry::hashing::Real6 loop_transform, core::Size loop_size ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols;
	using namespace protocols::loophash;
	using namespace conformation;
	using namespace kinematics;
	using namespace numeric::geometry::hashing;
	using namespace optimization;
	using namespace id;

	TR << "Getting hashmap ... " << std::endl;
	LoopHashMap &hashmap = loop_hash_library->gethash( loop_size );
	std::vector < core::Size > leap_index_bucket;
	TR << "Radial lookup ... " << std::endl;
	hashmap.radial_lookup( 0, loop_transform, leap_index_bucket );
	core::Size example_index = leap_index_bucket[0];
	TR << "Get the actual strucure index (not just the bin index) << " << std::endl;

	LeapIndex cp = hashmap.get_peptide( example_index );
	// Retrieve the actual backbone structure
	BackboneSegment new_bs;
	loop_hash_library->backbone_database().get_backbone_segment( cp.index, cp.offset, hashmap.get_loop_size() , new_bs );

	return new_bs;
}

int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;
		using namespace protocols;
		using namespace protocols::loophash;
		using namespace conformation;
		using namespace kinematics;
		using namespace numeric::geometry::hashing;
		using namespace optimization;
		using namespace id;


		// initialize core
		devel::init(argc, argv);

		TR << "Testing..." << std::endl;
		core::pose::Pose sample_pose;

		std::string seq = "TAKESMEFCKANDSMSHITMAKEAFCKSHITSTACKAFCKSHITSTACK";
		make_pose_from_sequence(sample_pose, seq , core::chemical::FA_STANDARD );

		// make an alpha helix
		for ( core::Size ir = 1; ir < sample_pose.size(); ir ++ ) {
			sample_pose.set_phi( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
			sample_pose.set_psi( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
			sample_pose.set_omega( ir, numeric::random::rg().uniform()*360.0 - 180.0 );
		}

		sample_pose.dump_pdb("testout.pdb");

		TR << "Importing test pose ..." << std::endl;
		// if(  option[ in::file::s ]().size() != 1){
		//  TR.Error << "Please supply a single inputfile with -in:file:s " << std::endl;
		// }
		// core::import_pose::pose_from_file( sample_pose, option[ in::file::s ]()[1] , core::import_pose::PDB_file);

		utility::vector1 < core::Size > loop_sizes;
		loop_sizes.push_back( 10 );
		loophash::LoopHashLibraryOP loop_hash_library = new loophash::LoopHashLibrary( loop_sizes, 0 , 1 );

		loop_hash_library->test_saving_library( sample_pose, 14, true /*deposit structure*/ );

		TR << "Save Library: " << std::endl;

		loop_hash_library->save_db();

		TR << "Reloading Library: " << std::endl;

		loophash::LoopHashLibraryOP read_loop_hash_library = new loophash::LoopHashLibrary( loop_sizes, 0 , 1 );
		read_loop_hash_library->load_db();

		TR << "Retesting Library: " << std::endl;

		read_loop_hash_library->test_saving_library( sample_pose, 14, false /*dont deposit structure*/ ) ;

		TR << "Done test" << std::endl;
		return -1;

		// add pose to the backbone database

		// Ok, take an example loop

		BackboneSegment pose_bs;
		core::Size ir = 16;
		core::Size loop_size = 15;
		TR << "Taking a sample loop from pose: " << ir <<  "  " << loop_size << std::endl;
		pose_bs.read_from_pose( sample_pose, ir, loop_size );
		Real6 loop_transform;
		if ( !get_rt_over_leap( sample_pose, ir, ir + loop_size, loop_transform ) ) {
			std::cerr << "Something went wrong" << std::endl;
			return -1;
		}


		BackboneSegment new_bs = get_backbone_segment( loop_hash_library,  loop_transform, loop_size );
		new_bs.print();

		loophash::LoopHashLibraryOP loop_hash_library_loaded = new loophash::LoopHashLibrary( loop_sizes );
		loop_hash_library_loaded->load_db();

		BackboneSegment new_bs_loaded = get_backbone_segment( loop_hash_library,  loop_transform, loop_size );
		new_bs_loaded.print();


	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


