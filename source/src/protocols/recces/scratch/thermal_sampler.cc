// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/recces/sampler/rna/MC_RNA_Suite.hh>
#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>
#include <protocols/recces/util.hh>
#include <protocols/recces/scratch/thermal_sampler.hh>

#include <core/id/TorsionID.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <utility/excn/Exceptions.hh>


// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>
#include <basic/options/keys/recces.OptionKeys.gen.hh>


using namespace core::pose;
using namespace basic::options;


using namespace core;
using namespace protocols;
using namespace protocols::recces;
using namespace protocols::moves;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::recces;
using namespace protocols::recces;
using utility::vector1;

namespace protocols {
namespace recces {
namespace scratch {


//////////////////////////////////////////////////////////////////////////////
utility::vector1<core::Real> get_torsions(
	utility::vector1<core::id::TorsionID> & torsion_ids,
	const Pose & pose
) {
	utility::vector1<core::Real> curr_torsions;
	for ( Size i = 1; i <= torsion_ids.size(); ++i ) {
		curr_torsions.push_back( pose.torsion( torsion_ids[i] ) );
	}
	return curr_torsions;
}

//////////////////////////////////////////////////////////////////////////////
void set_gaussian_stdevs(
	utility::vector1<protocols::recces::sampler::rna::MC_RNA_KIC_SamplerOP> & internal_bb_sampler,
	utility::vector1<protocols::recces::sampler::MC_OneTorsionOP> & chi_sampler,
	sampler::rna::MC_RNA_MultiSuite & standard_bb_sampler,
	moves::SimulatedTempering const & tempering,
	Size const & total_rsd,
	Size const & sampled_rsd,
	utility::vector1<bool> is_free
) {
	Real const temp( tempering.temperature() );
	Real internal_bb_stdev( 0.1 * pow( temp, 0.25 ) + 0.1);
	Real free_chi_stdev( 55 * pow( temp, 0.5 ) + 50 );
	Real chi_stdev( 5 * pow( temp , 0.5) + 15 );
	Real standard_bb_stdev( 8 * pow( temp, 0.5 ) / ( 2 * total_rsd + sampled_rsd ) );
	if ( temp < 0 ) {
		internal_bb_stdev = 0.5 ;
		free_chi_stdev = -1 ;
		chi_stdev = -1 ;
		standard_bb_stdev = -1 ;
	}
	for ( Size i = 1; i <= internal_bb_sampler.size(); ++i ) {
		internal_bb_sampler[i]->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( Size i = 1; i <= chi_sampler.size(); ++i ) {
		if ( is_free[i] ) {
			chi_sampler[i]->set_gaussian_stdev( free_chi_stdev );
		} else chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
	standard_bb_sampler.set_gaussian_stdev( standard_bb_stdev );
}

} //scratch
} //recces
} //protocols

