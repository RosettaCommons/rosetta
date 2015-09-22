// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopHashSampler.cc
/// @brief
/// @author Mike Tyka

#include <protocols/loophash/LoopHashRelaxProtocol.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/string_util.hh>
#include <protocols/relax/FastRelax.hh>
//#include <protocols/match/Hit.fwd.hh>
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
#include <protocols/loops/Loops.hh>

#include <numeric/geometry/hashing/SixDHasher.hh>
#include <numeric/random/random.hh>

// C++ headers
//#include <cstdlib>

#include <iostream>
#include <string>
#include <cstdio>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/vector1.hh>


// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>
#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


static THREAD_LOCAL basic::Tracer TR( "LocalHashRelaxProtocol" );


namespace protocols {
namespace loophash {

LoopHashRelaxProtocol::LoopHashRelaxProtocol( LoopHashLibraryOP library ) : library_(library)
{
	std::cout << "HERE!" << std::endl;
}

void
LoopHashRelaxProtocol::manual_call( core::pose::Pose& pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;

	core::Size skim_size = 2;

	if ( !library_ ) std::cout << std::string( "FATAL ERROR: Libary not loaded ") ;

	//core::Real mpi_metropolis_temp_ = option[ lh::mpi_metropolis_temp ]();

	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	LoopHashSampler  lsampler( library_, simple_inserter );
	lsampler.set_min_bbrms( option[ lh::min_bbrms ]()   );
	lsampler.set_max_bbrms( option[ lh::max_bbrms ]() );
	lsampler.set_min_rms( option[ lh::min_rms ]() );
	lsampler.set_max_rms( option[ lh::max_rms ]() );
	lsampler.set_max_struct( skim_size );

	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
	} else {
		native_pose = pose; // jsut make a copy of the current pose - rmses will be relative to starts
	}

	// Set up contraints
	ScoreFunctionOP fascorefxn = get_score_function();

	// convert pose to centroid pose:
	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);
	}


	// See if a loopfile was defined - if so restrict sampling to those loops
	//    protocols::loops::Loops loops(true);
	//    utility::vector1< core::Size > selection;
	//    loops.get_residues( selection );
	//    TR << "Userdefined Loopregions: " << loops.size() << std::endl;
	//    TR << loops << std::endl;
	//    TR << "Residues: ";
	//    for( core::Size i=1; i <= selection.size(); ++i) TR <<  selection[i] << " ";
	//    TR << std::endl;

	//read_coord_cst(); //include this function later !

	//core::Size total_starttime = time(NULL);

	static int casecount = 0;
	core::pose::Pose opose = pose;
	std::vector< core::io::silent::SilentStructOP > lib_structs;

	TR.Info << "Loophash apply function ! " << std::endl;

	//protocols::relax::FastRelax *qrelax = new protocols::relax::FastRelax( fascorefxn, 1 );
	protocols::relax::FastRelax *relax = new protocols::relax::FastRelax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() );

	// convert pose to centroid pose:
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID);
	core::pose::set_ss_from_phipsi( pose );

	// Generate alternate structures
	core::Size starttime2 = time(NULL);
	core::Size sampler_chunk_size = 1;
	core::Size start_res;
	core::Size stop_res;
	do {
		start_res = std::max( core::Size(2), core::Size(rand()%(pose.total_residue() - sampler_chunk_size - 2 )) );
		stop_res  = std::min( core::Size(pose.total_residue()), core::Size(start_res + sampler_chunk_size - 1 )  );

		// If a loopfile was given choose your insertion site from there
		//    TR.Info << "Selection size: " << selection.size() << std::endl;
		//    if( selection.size() > 0 ){
		//      utility::vector1< core::Size > temp_selection = selection;
		//      std::random_shuffle( temp_selection.begin(), temp_selection.end());
		//      start_res = std::max( core::Size(2), core::Size( temp_selection[1] ) );
		//      stop_res  = std::min( core::Size(pose.total_residue()), core::Size(start_res + sampler_chunk_size - 1)  );
		//      TR.Info << "SubselectionSample: " << start_res << " - " << stop_res << std::endl;
		//    }

		lsampler.set_start_res( start_res );
		lsampler.set_stop_res(  stop_res );
		lsampler.build_structures( pose, lib_structs );
		core::Size endtime2 = time(NULL);
		//core::Size loophash_time = endtime2 - starttime2;
		TR.Info << "FOUND (" << start_res << " to " << stop_res << "): "
			<< lib_structs.size() << " states in time: "
			<< endtime2 - starttime2 << " s " << std::endl;

		// try again if we have failed to find structures << danger this is totally an infinte loop here

	} while ( lib_structs.size() == 0 );

	// choose up to "skim_size" of them
	std::random_shuffle( lib_structs.begin(), lib_structs.end());
	std::vector< core::io::silent::SilentStructOP > select_lib_structs;
	for ( core::Size k=0; k< std::min<core::Size>(skim_size, lib_structs.size() ) ; k++ ) {
		select_lib_structs.push_back( lib_structs[k] );
	}

	core::Real bestcenscore = MAXIMAL_FLOAT;
	//core::Size bestcenindex = 0;
	for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
		core::pose::Pose rpose;
		select_lib_structs[h]->fill_pose( rpose );

		//rpose.dump_pdb("struct_" + string_of(h) + ".pdb" );

		core::Real refrms = 0;
		core::Real rms_factor = 10.0;
		core::Real decoy_score = select_lib_structs[h]->get_energy("lh_censcore") + refrms * rms_factor;

		select_lib_structs[h]->add_energy( "refrms",     refrms,      1.0 );
		select_lib_structs[h]->add_energy( "comb_score", decoy_score, 1.0 );
		TR.Info << "refrms: " << refrms << "  Energy: " << decoy_score << std::endl;
		if ( decoy_score < bestcenscore ) {
			bestcenscore = decoy_score;
			//bestcenindex = h;  // set but never used ~Labonte
		}
	}
	TR.Info << "Best:" << "  Energy: " << bestcenscore << std::endl;

	/// For fullatom goodness, continue
	core::Size starttime = time(NULL);
	relax->batch_apply( select_lib_structs );
	core::Size endtime = time(NULL);
	//core::Size batchrelax_time = endtime - starttime;
	TR.Info << "Batchrelax time: " << endtime - starttime << " for " << select_lib_structs.size() << " structures " << std::endl;

	core::Real bestscore = MAXIMAL_FLOAT;
	//core::Size bestindex = 0;
	core::pose::Pose relax_winner;
	for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
		TR.Info << "DOING: " << h << " / " << select_lib_structs.size() << std::endl;
		core::pose::Pose rpose;

		select_lib_structs[h]->fill_pose( rpose );

		//core::Real score = CA_rmsd( native_pose, rpose );
		core::Real score = (*fascorefxn)(rpose);
		TR.Info << "score: " << h << "  " << score << std::endl;

		select_lib_structs[h]->add_energy("lh_score_new", score );

		if ( score < bestscore ) {
			bestscore = score;
			//bestindex = h;  // set but never used ~Labonte
			relax_winner = rpose;
		}
	}
	casecount++;


	std::cout << "SETTING Relax winner" << std::endl;
	pose = relax_winner;

}


void
LoopHashRelaxProtocol::apply( core::pose::Pose& /*pose*/ )
{
}

protocols::moves::MoverOP LoopHashRelaxProtocol::clone() const {
	return protocols::moves::MoverOP( new LoopHashRelaxProtocol( *this ) );
}

protocols::moves::MoverOP LoopHashRelaxProtocol::fresh_instance() const {
	return protocols::moves::MoverOP( new LoopHashRelaxProtocol( library_ ) );
}

} // namespace loops
} // namespace protocols
