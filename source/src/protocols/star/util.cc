// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <string>

// External headers
#include <boost/format.hpp>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/util/SwitchResidueTypeSet.hh>

namespace protocols {
namespace star {

void emit_intermediate(const core::pose::Pose& pose, const std::string& silent_filename) {
	using core::Size;
	using core::io::silent::SilentFileData;
	using core::io::silent::SilentStructFactory;
	using core::io::silent::SilentStructOP;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[OptionKeys::abinitio::debug]() ) {
		return;
	}

	static Size num_structures = 0;

	SilentStructOP silent = SilentStructFactory::get_instance()->get_silent_struct_out();
	silent->fill_struct(pose, str(boost::format("model_%d") % ++num_structures));

	SilentFileData sfd;
	sfd.write_silent_struct(*silent, silent_filename);
}

void simple_fold_tree(core::pose::Pose & pose) {
	pose.fold_tree(core::kinematics::FoldTree(pose.total_residue()));
}

void to_centroid(core::pose::Pose & pose) {
	if ( !pose.is_centroid() ) {
		core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
	}
}

}  // namespace star
}  // namespace protocols
