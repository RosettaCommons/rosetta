// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/util.hh> // for reroot
#include <protocols/stepwise/modeler/rna/util.hh> // for virtualize_free_rna_moieties
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh> // for fill_vdw_cached_rep_screen_info_from_command_line
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <utility/stream_util.hh>
#include <utility/vector1.functions.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.setup.FullModelInfoSetupFromCommandLine" );

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;
using namespace core::id;
using utility::vector1;
using std::pair;
using std::make_pair;
using core::kinematics::FoldTree;

//////////////////////////////////////////////////////////////////////////////////////
//
// All setup functions, including PDB readin & FullModelInfo initialization for
//  stepwise application.
//
// Getting to be a big file, might be good to separate some functions out, and/or
//  turn into a class?
//
//                                -- rhiju, 2014
//
//  TODO: Make a class or mover?
//   -- AMW: few of these functions have shared data that merit class-hood, but
//      I agree that they should be grouped in a manner more sane than just
//      grouping into a file of related functions!
//  TODO: Move outside protocols/stepwise (since will soon be in use in fragment assembly)
//   -- AMW: most of these functions are now in pose/util!
//  TODO: Set up -cyclize flag (like what it is available for RNA in farna)
//   -- done!
//
//////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace setup {

// AMW TODO: version that takes an OptionsCollection
// AMW: for now, create version that leaves
// alone if it has already been initialized.
///////////////////////////////////////////////////////////////
void
initialize_native_and_align_pose( PoseOP & native_pose,
	PoseOP & align_pose,
	core::chemical::ResidueTypeSetCAP rsd_set,
	PoseCOP start_pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( ( !native_pose || native_pose->size() == 0 ) && option[ in::file::native ].user() )  {
		// If native is on the command line, stuff all that in align_pose.
		align_pose = native_pose = core::import_pose::get_pdb_with_full_model_info( option[ in::file::native ](), rsd_set );
	} else if ( ( !align_pose || align_pose->size() == 0 ) && native_pose ) {
		// If native wasn't specified on the command line but it is already set up
		// nonetheless...
		align_pose = native_pose;
	}

	if ( option[ OptionKeys::stepwise::new_align_pdb ].user() ) {
		align_pose = core::import_pose::get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::new_align_pdb ](), rsd_set );
	} else if ( ( !align_pose || align_pose->size() == 0 )  && option[ OptionKeys::stepwise::align_pdb ].user() ) {
		align_pose = core::import_pose::get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::align_pdb ](), rsd_set );
	}

	if ( align_pose == 0 && option[ in::file::s ].user() ) {
		align_pose = start_pose->clone();
	}
	if ( option[ OptionKeys::stepwise::virtualize_free_moieties_in_native ]() ) { // could generalize to proteins
		if ( native_pose != nullptr )  modeler::rna::virtualize_free_rna_moieties( *native_pose );
		if ( align_pose  != nullptr ) modeler::rna::virtualize_free_rna_moieties( *align_pose );
	}
}


} //setup
} //stepwise
} //protocols
