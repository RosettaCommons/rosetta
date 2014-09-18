// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/FastGapMover.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/loophash/FastGapMover.hh>

#include <basic/Tracer.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/cm.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>

#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>


#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// C++ headems
// AUTO-REMOVED #include <ctime>
#include <algorithm>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>



namespace protocols {
namespace loophash {


using namespace core;
using namespace io;
using namespace pdb;
using namespace chemical;
using namespace conformation;
using namespace protocols::loophash;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace utility;

using core::pose::PoseOP;
using core::pose::PoseAP;
using core::pose::Pose;
using utility::vector1;
using core::Size;
using std::string;

static thread_local basic::Tracer TR( "protocols.loophash.FastGapMover" );

FastGapMover::FastGapMover() :
		Mover("FastGapMover"),
		min_loop_size_(3),
		max_loop_size_(14),
		max_rms_(0.3),
		min_rms_(0.0),
		non_ideal_(true) 
	{
	// initialize lhlibrary
	utility::vector1 < core::Size > loop_sizes;
	for (Size i = min_loop_size_; i<= max_loop_size_; i++ ) {
		loop_sizes.push_back(i);
	}

	simple_inserter_ = new LocalInserter_SimpleMin();
	lhlibrary_ = new LoopHashLibrary( loop_sizes );
	lhlibrary_->load_mergeddb();

	lhsampler_ = new LoopHashSampler( lhlibrary_, simple_inserter_ );
} 

void
FastGapMover::apply( pose::Pose & pose ) {
		lhsampler_->set_max_rms( max_rms_ );
		lhsampler_->set_nonideal( non_ideal_ );
		lhsampler_->set_min_rms( min_rms_ );

		// copy pose
		PoseOP working_pose = new Pose(pose);
		//Size idx = 1;

		// convert pose to centroid pose:
		if( working_pose->is_fullatom() ){
			core::util::switch_to_residue_type_set( *working_pose, core::chemical::CENTROID);
		}
		// Now go through each gap and try increasingly larger lh until something is returned
		Size next_gap = 0;
		Real gap_dist;
		find_next_gap( *working_pose, next_gap, gap_dist );
		while( next_gap != 0 ) {
			TR << "Attempting to fix gap following residue " << next_gap << std::endl;
			std::vector< Pose > lib_structs;
			
			// gogo loophash
			// increase loophash size until we get anything returned
			//Size loop_size = std::max((Size)(gap_dist/3.5), min_loop_size_); // no point in trying anything that cant reach across the gap
			Size loop_size = min_loop_size_; // no point in trying anything that cant reach across the gap
			while( lib_structs.size() == 0 && loop_size < max_loop_size_ ) {
				TR << "Trying loopsize " << ++loop_size << std::endl;
				lhsampler_->set_start_res( next_gap + 3 < loop_size ? 0 : next_gap + 3 - loop_size );
				lhsampler_->set_stop_res ( next_gap );
				lhsampler_->close_gaps( *working_pose, lib_structs, loop_size );
			}
			if( lib_structs.size() != 0 ) {
				working_pose = new Pose (lib_structs[0]);
			}

			find_next_gap( *working_pose, next_gap, gap_dist );
		}
		pose = *working_pose;

} // apply

// lifted straight from protocols/idealize/idealize.cc
void
FastGapMover::find_next_gap( Pose & pose, Size & idx, Real & gap_distance ) {
	// squared distance at which bond is considered discontinuous
	Real const chain_break_cutoff = { 4.0 };
	Size const nres ( pose.total_residue() );

	// find chain breaks to add to gaplist
	kinematics::FoldTree f( pose.fold_tree() );
	for ( Size i = idx + 1; i < nres; ++i ) {
			//bool chain_break = false;
			Size j = i+1;
			conformation::Residue const & rsd = pose.residue(i);
			conformation::Residue const & next_rsd = pose.residue(j);
			if (rsd.is_polymer() && next_rsd.is_polymer()) {
					Real dist_squared = rsd.atom( rsd.upper_connect_atom() ).xyz().distance_squared(next_rsd.atom( next_rsd.lower_connect_atom() ).xyz());
					gap_distance = std::sqrt(dist_squared);
					if (dist_squared > chain_break_cutoff || dist_squared < 0.1) {
							//chain_break = true;  // set but never used ~Labonte
							idx = i;
							return;
					}
			}
	}
	idx = 0;
}
string
FastGapMover::get_name() const {
	return "FastGapMover";
}


} // namespace idealize
} // namespace protocols
