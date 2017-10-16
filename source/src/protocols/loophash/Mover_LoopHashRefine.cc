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

// Must place this on top because of an ambiguous symbol issue with byte
#ifdef BOINC
	#include <protocols/boinc/boinc.hh>  // REQUIRED FOR WINDOWS
#endif

// libRosetta headers
#include <protocols/loophash/Mover_LoopHashRefine.hh>


#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/relax/FastRelax.hh>
// #include <protocols/match/Hit.fwd.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <protocols/moves/Mover.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/optimization/MinimizerOptions.hh>

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <numeric/random/random.hh>

// C++ headers
//#include <cstdlib>

#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>


#ifdef WIN32
#include <ctime>
#endif


static THREAD_LOCAL basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace numeric::geometry::hashing;
using namespace protocols::frag_picker;
using namespace protocols::loophash;


namespace protocols {
namespace loophash {


void
Mover_LoopHashRefine::apply( core::pose::Pose& pose )
{
	if ( !library_ ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//std::string prefix = option[ out::prefix ]();
	core::Size skim_size = option[ lh::skim_size ]();

	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	LoopHashSampler  lsampler( library_, simple_inserter );


	ScoreFunctionOP fascorefxn = core::scoring::get_score_function();
	protocols::relax::FastRelax relax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() );
	core::pose::PoseOP native_pose;
	if (  option[ in::file::native ].user() ) {
		native_pose = core::pose::PoseOP( new Pose );
		core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		relax.set_native_pose( native_pose );
	}

	for ( int round = 0; round < option[ OptionKeys::lh::rounds ]; round ++ ) {
		std::string checkpoint_id = "chk" + ObjexxFCL::string_of( round );
		if ( !checkpoints_.recover_checkpoint( pose, get_current_tag(), checkpoint_id, true, true ) ) {
			core::pose::Pose opose = pose;
			std::vector< core::io::silent::SilentStructOP > lib_structs;

			// convert pose to centroid pose and apply loophasher
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
			core::pose::set_ss_from_phipsi( pose );
			core::Size starttime2 = time(nullptr);

			core::Size lcount = 0;
			while ( lib_structs.size() < skim_size ) {
				core::Size resnum = numeric::random::random_range(1,pose.size()-2);

				lsampler.set_start_res ( resnum );
				lsampler.set_stop_res ( resnum );
				lsampler.build_structures( pose, lib_structs );
				lcount++;
				TR.Info << "Lcount: " << lcount << std::endl;
			}
			core::Size endtime2 = time(nullptr);
			TR.Info << "FOUND " << lib_structs.size() << " alternative states in time: " << endtime2 - starttime2 << std::endl;

			// write out the centroid structures if desired
			if ( (  option[ OptionKeys::lh::write_centroid_structs ]() ) ||
					(  option[ OptionKeys::lh::centroid_only ]() ) ) {
				core::io::silent::SilentFileOptions opts;
				core::io::silent::SilentFileData sfd( opts );
				std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
				silent_file_ += ".centroid.out" ;
				for ( core::Size h = 0; h < lib_structs.size(); h++ ) {
					core::pose::Pose rpose;
					lib_structs[h]->fill_pose( rpose );
					lib_structs[h]->add_energy( "round", round, 1.0 );
					if ( native_pose ) {
						core::Real rms = scoring::CA_rmsd( *native_pose, rpose );
						core::Real gdtmm = scoring::CA_gdtmm( *native_pose, rpose );
						lib_structs[h]->add_energy( "rms", rms, 1.0 );
						lib_structs[h]->add_energy( "gdtmm", gdtmm, 1.0 );
					}
					lib_structs[h]->set_decoy_tag( "S_" + ObjexxFCL::string_of( round ) + "_" + ObjexxFCL::string_of(  h )  );
					lib_structs[h]->sort_silent_scores();
					sfd.write_silent_struct( *(lib_structs[h]) , silent_file_ );
				}
			}

			// abort if centroid is all we're doing
			if ( option[ OptionKeys::lh::centroid_only ]() ) break;

			// Choose a set of structures to refine/relax
			//std::random__shuffle( lib_structs.begin(), lib_structs.end());
			numeric::random::random_permutation(lib_structs.begin(), lib_structs.end(), numeric::random::rg());

			std::vector< core::io::silent::SilentStructOP > select_lib_structs;
			for ( core::Size k=0; k< std::min<core::Size>(skim_size, lib_structs.size() ) ; k++ ) {
				select_lib_structs.push_back( lib_structs[k] );
			}


			// Batch relax the result:
			core::Real bestscore = MAXIMAL_FLOAT;
			core::Size bestindex = 0;
			core::Size starttime = time(nullptr);
			relax.batch_apply( select_lib_structs );
			core::Size endtime = time(nullptr);
			TR.Info << "Batchrelax time: " << endtime - starttime << " for " << select_lib_structs.size() << " structures " << std::endl;


			for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
				TR.Info << "DOING: " << h << " / " << select_lib_structs.size() << std::endl;
				core::pose::Pose rpose;
				select_lib_structs[h]->fill_pose( rpose );
				core::Real score = (*fascorefxn)(rpose);
				TR.Info << "score: " << h << "  " << score << std::endl;
				if ( score < bestscore ) {
					bestscore = score;
					bestindex = h;
					pose = rpose;
				}
			}

			core::io::silent::SilentFileOptions opts;
			core::io::silent::SilentFileData sfd( opts );
			std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
			for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
				if ( h == bestindex ) {
					core::pose::Pose rpose;
					select_lib_structs[h]->fill_pose( rpose );
#ifdef BOINC_GRAPHICS
					//mjo moving 'score' into ifdef to remove unused variable warning
					core::Real score = select_lib_structs[h]->get_energy("score");
					boinc::Boinc::update_graphics_low_energy( pose, score );
					boinc::Boinc::update_graphics_last_accepted( pose, score );
#endif
					select_lib_structs[h]->add_energy( "round", round, 1.0 );
					select_lib_structs[h]->set_decoy_tag( "S_" + ObjexxFCL::string_of( round ) + "_" + ObjexxFCL::string_of(  h )  );
					select_lib_structs[h]->sort_silent_scores();

					sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
				}
			}
			checkpoints_.checkpoint( pose, get_current_tag(), checkpoint_id,  true );
		}

		// if in boinc mode let boinc break out prematurely! This will make user's runtimes near perfect and users happy
#ifdef BOINC
			if (protocols::boinc::Boinc::worker_is_finished( round + 1 )){
				std::cerr << "BOINC Finishing normally." << std::endl;
				break; // we're done no matter what nstruct syays
			}
#endif
	}

	TR << "Finished serial loophash." << std::endl;

}


// this is in protocols so that it can be caled from the boinc main function
int loophash_main(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::loophash;

	utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
	LoopHashLibraryOP loop_hash_library( new LoopHashLibrary( loop_sizes ) );

	// Run simple sampling run test or create the db ?
	if ( option[lh::create_db]() ) { ;
		loop_hash_library->create_db();
		loop_hash_library->save_db();
		return 0;
	}

	Mover_LoopHashRefineOP lh_sampler( new Mover_LoopHashRefine( loop_hash_library ) );

	// Normal mode with external loophash library
	loop_hash_library->load_db();
	try{
		protocols::jd2::JobDistributor::get_instance()->go( lh_sampler );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
	}

	return 0;
}


}
}
