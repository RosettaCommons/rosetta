// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/wum/WorkUnitBase.cc
/// @brief
/// @author Mike Tyka

#include <protocols/loophash/WorkUnit_LoopHash.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/LoopHashSampler.hh>

//Auto Headers
#include <utility/vector1.hh>
//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

namespace protocols {
namespace loophash {


static THREAD_LOCAL basic::Tracer TR( "WorkUnit_LoopHash" );


WorkUnit_LoopHash::WorkUnit_LoopHash( core::Size start_ir, core::Size end_ir, core::Size ssid ):
	WorkUnit_SilentStructStore()
{
	set_defaults();
	set_start(start_ir);
	set_end(end_ir);
	set_ssid(ssid);
}

void
WorkUnit_LoopHash::set_defaults()
{
	TR.Debug << "Setting type to WU_Type_LoopHash" << std::endl;
}

void
WorkUnit_LoopHash::init_from_cmd( const core::Size mpi_rank )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1 < core::Size > loop_sizes = option[ OptionKeys::lh::loopsizes]();
	core::Size num_partitions = option[ OptionKeys::wum::n_slaves_per_master]();
	if ( option[ OptionKeys::lh::num_partitions].user() ) {
		num_partitions = option[ OptionKeys::lh::num_partitions]();
	}
	core::Size assigned_num = mpi_rank % num_partitions;
	try{
		library_ = LoopHashLibraryOP( new LoopHashLibrary( loop_sizes, num_partitions, assigned_num ) );
		// load initial library from disk
		library_->load_mergeddb();
	}
catch( utility::excn::EXCN_Msg_Exception & e ){
	e.show( std::cout );
	e.show( std::cerr );
	throw;
}
	library_->mem_foot_print();
}

void
WorkUnit_LoopHash::run()
{
	using namespace core::pose;
	using namespace protocols::loops;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	core::io::silent::SilentStructCOP start_struct = decoys().get_struct(0);
	core::pose::Pose pose;
	decoys().get_pose( 0, pose );

	// clear the sotre of structures
	decoys().clear();

	runtime_assert( library_ != nullptr );

	TR << "Executing WorkUnit_LoopHash_Mover..." << std::endl;
	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	LoopHashSampler  lsampler( library_, simple_inserter );
	lsampler.set_start_res( get_start()  );
	lsampler.set_stop_res ( get_end() );
	lsampler.set_min_bbrms( 20.0   );
	lsampler.set_max_bbrms( 1400.0 );
	lsampler.set_min_rms( 0.5 );
	lsampler.set_max_rms( 4.0 );
	if ( option[ OptionKeys::lh::bss]() ) lsampler.set_nonideal(true);

	// convert pose to centroid pose:
	if ( pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
	}
	core::pose::set_ss_from_phipsi( pose );

	TR.Info << "Running loophash function: Start: " << get_start() << " End: " << get_end() << std::endl;
	core::Size starttime = time(nullptr);

	lsampler.build_structures( pose, decoys().store() );

	core::Size endtime = time(nullptr);
	TR.Info << "Build " << decoys().size() << " structures in " << endtime - starttime << " s " << std::endl;

	// transfer any tags from input structure:

	const core::io::silent::SilentStruct *ss2 = &( *start_struct );

	for ( auto & it : decoys() ) {

		core::io::silent::SilentStruct *ss = &(*it);
		// preserve the centroid score!
		core::Real censcore = it->get_energy("censcore");
		it->copy_scores( *start_struct );
		it->add_energy("censcore", censcore );
		std::string new_usid = protocols::wum::generate_unique_structure_id();
		it->add_string_value( "husid", it->get_string_value("husid") + "." + new_usid );
		it->add_string_value( "usid", new_usid );
		it->add_energy( "state", 1 );

		if ( option[ OptionKeys::lh::bss]() ) {
			const core::io::silent::BinarySilentStruct *pss2 = dynamic_cast< const core::io::silent::BinarySilentStruct* > ( ss2 );
			core::io::silent::BinarySilentStruct *pss = dynamic_cast< core::io::silent::BinarySilentStruct* > ( ss );
			if ( (pss != nullptr) && (pss2 != nullptr) ) {
				TR.Debug << "LoophashResult: " << pss->CA_rmsd( *pss2 ) << "  " <<  pss->get_energy("censcore") << std::endl;
			} else {
				TR << "LoophashResult: ERROR, dynamic cast to BSS failed"  << std::endl;
			}
		} else {
			const core::io::silent::ProteinSilentStruct *pss2 = dynamic_cast< const core::io::silent::ProteinSilentStruct* > ( ss2 );
			core::io::silent::ProteinSilentStruct *pss = dynamic_cast< core::io::silent::ProteinSilentStruct* > ( ss );
			if ( (pss != nullptr) && (pss2 != nullptr) ) {
				TR.Debug << "LoophashResult: " << pss->CA_rmsd( *pss2 ) << "  " <<  pss->get_energy("censcore") << std::endl;
			} else {
				TR << "LoophashResult: ERROR, dynamic cast to PSS failed"  << std::endl;
			}
		}
	}
}


}
}

