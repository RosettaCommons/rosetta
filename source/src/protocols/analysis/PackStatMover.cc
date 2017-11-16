// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   PackStatMover.cc
///
/// @brief
/// @author Monica Berrondo

#include <protocols/analysis/PackStatMover.hh>

#include <core/scoring/packstat/compute_sasa.hh>

//#include <protocols/jobdist/standard_mains.hh>

#include <core/pose/PDBInfo.hh>

#include <basic/options/option.hh>


#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/packstat.OptionKeys.gen.hh>

#include <utility/vector1.hh>


//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.analysis.packstat" );
using namespace core;

namespace protocols {
namespace analysis {

PackStatMover::PackStatMover() : moves::Mover() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// bool verbose_ = verbose;
	packstat_pdb_                     = option[ OptionKeys::packstat::packstat_pdb ]();
	include_water_                    = option[ OptionKeys::packstat::include_water ]();
	surface_accessibility_            = option[ OptionKeys::packstat::surface_accessibility ]();
	residue_scores_                   = option[ OptionKeys::packstat::residue_scores ]();
	cavity_burial_probe_radius_       = option[ OptionKeys::packstat::cavity_burial_probe_radius ]();
	oversample_                       = option[ OptionKeys::packstat::oversample ]();

}


void
PackStatMover::apply( pose::Pose & pose )
{
	using namespace core::scoring::packstat;
	using namespace std;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL::format;
	using namespace numeric;
	using namespace utility;

	std::string fname = pose.pdb_info()->name();


	PosePackData pd = pose_to_pack_data( pose, include_water_ );

	core::Real packing_score = compute_packing_score( pd, oversample_ );

	TR << "packing score: " << fname << " " << packing_score;
	if ( include_water_ ) {
		TR << " ( with waters? (pose?) )";
	}
	TR << std::endl;

	utility::vector1<core::Real> res_scores; // needed if output res scores or pdb
	if ( packstat_pdb_ || residue_scores_ ) {
		res_scores = compute_residue_packing_scores( pd, oversample_ );
	}

	// if( packstat_pdb_ ) {
	//  output_packstat_pdb( fname, res_scores );
	// }

	if ( residue_scores_ ) {
		for ( int i = 1; i <= (int)res_scores.size(); ++i ) {
			TR << "packing score: residue " << i << " " << res_scores[i] << std::endl;
		}
	}


}

std::string
PackStatMover::get_name() const {
	return "PackStatMover";
}

} // analysis
} // protocols
