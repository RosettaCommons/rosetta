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
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>

#include <protocols/recces/sampler/rna/MC_RNA_MultiSuite.hh>
#include <protocols/moves/SimulatedTempering.hh>
#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/recces/scratch/thermal_sampler.hh>

#include <core/id/TorsionID.fwd.hh>
#include <protocols/recces/sampler/MC_OneTorsion.hh>

// C++ headers

// option key includes



// option key includes

#include <basic/options/keys/recces.OptionKeys.gen.hh>

#include <cmath> // MANUAL IWYU

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
	for ( core::Size i = 1; i <= torsion_ids.size(); ++i ) {
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
	core::Size const & total_rsd,
	core::Size const & sampled_rsd,
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
	for ( core::Size i = 1; i <= internal_bb_sampler.size(); ++i ) {
		internal_bb_sampler[i]->set_gaussian_stdev( internal_bb_stdev );
	}
	for ( core::Size i = 1; i <= chi_sampler.size(); ++i ) {
		if ( is_free[i] ) {
			chi_sampler[i]->set_gaussian_stdev( free_chi_stdev );
		} else chi_sampler[i]->set_gaussian_stdev( chi_stdev );
	}
	standard_bb_sampler.set_gaussian_stdev( standard_bb_stdev );
}

} //scratch
} //recces
} //protocols

