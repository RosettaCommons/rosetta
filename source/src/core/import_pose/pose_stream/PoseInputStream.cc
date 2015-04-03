// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseInputStream.cc
/// @brief
/// @author James Thompson

#include <core/import_pose/pose_stream/PoseInputStream.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/io/raw_data/DisulfideFile.hh>

#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

// STL
#include <utility>

namespace core {
namespace import_pose {
namespace pose_stream {

void PoseInputStream::preprocess_pose( core::pose::Pose & pose ) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	if ( option[ james::debug ]() ) return;

	core::pose::initialize_disulfide_bonds(pose);

	// add constraints if specified by the user.
	// do this in a mover instead!
	//using namespace core::scoring::constraints;
	//if ( option[ constraints::cst_file ].user() ) {
	//	core::scoring::constraints::ConstraintSetOP
	//	cstset_ = ConstraintIO::get_instance()->read_constraints(
	//		get_cst_file_option(), new ConstraintSet, pose
	//	);
	//	pose.constraint_set( cstset_ );
	//}
} // PoseInputStream::preprocess_pose

utility::vector1< core::pose::PoseOP > PoseInputStream::get_all_poses(
	core::chemical::ResidueTypeSet const & residue_set
) {
	utility::vector1< core::pose::PoseOP > pose_list;
	reset();
	while ( has_another_pose() ) {
		core::pose::PoseOP pose( new core::pose::Pose );
		fill_pose( *pose, residue_set );
		pose_list.push_back( pose );
	}
	return pose_list;
}

} // pose_stream
} // import_pose
} // core
