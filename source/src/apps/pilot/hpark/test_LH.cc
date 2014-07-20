// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Project
#include <apps/pilot/hpark/NMsearch.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

// Pose
#include <core/id/AtomID.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// LoopHash
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <core/io/silent/SilentStruct.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// MD
#include <protocols/md/CartesianMD.hh>

// etc.
#include <numeric/random/random.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <sys/time.h>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::loophash;

static numeric::random::RandomGenerator RG( 12322 ); //Magic number??

protocols::loophash::LoopHashSamplerOP
setup_LHsampler()
{
  utility::vector1 < core::Size > loop_sizes;
  loop_sizes.push_back( 10 );

  LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes, 0, 1 );
  loop_hash_library->load_db();

	LocalInserter_CartesianOP simple_inserter( new LocalInserter_Cartesian() );
	protocols::loophash::LoopHashSamplerOP LHsampler = new LoopHashSampler::LoopHashSampler( loop_hash_library, simple_inserter );

	LHsampler->set_conserve_frame( true );
  LHsampler->set_min_bbrms( 20.0 );
  LHsampler->set_max_bbrms( 1400.0 );
  LHsampler->set_min_rms( 0.5 );
  LHsampler->set_max_rms( 4.0 );
	LHsampler->set_nonideal( true );

	LHsampler->set_max_struct( 1 );
	return LHsampler;
}

void
loophash_search( pose::Pose const pose,
								 LoopHashSamplerOP &LHsampler,
								 Size const ires )
{

	chemical::ResidueTypeSetCAP cen_rsdset =
		chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	protocols::moves::MoverOP tocen
		= new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );

	// Pick replacing resid
	LHsampler->set_start_res( ires );
	LHsampler->set_stop_res( ires );

	// Get centroid pose
	pose::Pose pose_cen( pose );
	tocen->apply( pose_cen );
	std::vector< core::io::silent::SilentStructOP > decoys;

	LHsampler->build_structures( pose_cen, decoys );

	/*
	// Get best among LoopHash results
	pose::Pose pose_work;
	pose::Pose pose_tmp( pose );
	for( Size i_decoy = 0; i_decoy < decoys.size(); ++i_decoy ){
		// Convert decoy to pose_work
		decoys[i_decoy]->fill_pose( pose_work, *cen_rsdset );
	}
	*/
}

int main( int argc, char *argv [] ){

  devel::init(argc, argv);

	using namespace myspace;

  core::pose::Pose pose, native;
  core::chemical::ResidueTypeSetCAP rsd_set
    = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

  core::import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ](1) );
  //core::import_pose::pose_from_pdb( native, *rsd_set, option[ in::file::native ]() );

	LoopHashSamplerOP LHsampler = setup_LHsampler();
	loophash_search( pose, LHsampler, 20 );

  return 0;
}

