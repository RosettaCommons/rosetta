// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/util.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/backrub/util.hh>

/////////// These were copied from the app.  Figure out what needs to go here before merging  //////////

// Protocols Headers


// Core Headers
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <platform/types.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <fstream>

static basic::Tracer TR("protocols.backrub.util");

namespace protocols {
namespace backrub {


bool
read_fold_tree_from_file(
	core::kinematics::FoldTree & foldtree,
	std::string filepath)
{
	std::ifstream filestream(filepath.c_str());

	while ( filestream.good() ) {

		std::string line;
		std::string key;

		getline(filestream, line);
		if ( filestream.fail() ) {
			//TR << "getline() failed" << std::endl;
			return false;
		}

		std::istringstream linestream(line);
		linestream >> key;
		if ( key == "FOLD_TREE" ) {
			linestream.clear();
			linestream.seekg(0, std::ios::beg);
			linestream >> foldtree;
			if ( linestream.fail() ) {
				TR << "FoldTree parsing failed" << std::endl;
				return false;
			} else {
				return true;
			}
		}
	}

	return false;
}

bool
read_fold_tree_from_file( core::pose::Pose & pose, std::string filepath)
{
	core::kinematics::FoldTree foldtree;

	if ( read_fold_tree_from_file(foldtree, filepath) ) {
		if ( foldtree.nres() == pose.size() ) {
			pose.fold_tree(foldtree);
			return true;
		} else {
			TR << "Different number of residues in Pose (" << pose.size() << ") and FoldTree (" << foldtree.nres()
				<< ")" << std::endl;
		}
	}

	return false;
}


void
append_fold_tree_to_file(
	core::kinematics::FoldTree const & foldtree,
	std::string file_path
)
{
	std::ofstream filestream( file_path.c_str(), std::ios::out|std::ios::app );
	if ( filestream.good() ) {
		filestream << foldtree << std::endl;
		filestream.close();
	} else {
		TR << "couldn't open file to append FoldTree" << std::endl;
	}
}


utility::vector1<core::Size>
positions_incompatible_with_task(
	core::pose::Pose & pose,
	core::pack::task::PackerTask & packertask)
{
	utility::vector1<core::Size> incompatible_positions;

	debug_assert(pose.size() == packertask.total_residue());

	// iterate over all residues to see if they're compatible
	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		// only check packable residues for compatibility
		if ( packertask.pack_residue(i) ) {

			// assume residue is incompatible
			bool incompatible(true);

			// check to see if pose residue type is in list of allowed residue types
			core::pack::task::ResidueLevelTask const & residueleveltask(packertask.residue_task(i));
			for ( auto iter(residueleveltask.allowed_residue_types_begin());
					iter != residueleveltask.allowed_residue_types_end(); ++iter ) {

				if ( (*iter)->name() == pose.residue_type(i).name() ) incompatible = false;
			}

			if ( incompatible ) incompatible_positions.push_back(i);
		}
	}

	return incompatible_positions;
}

utility::vector1<core::Size>
get_pivot_residues_from_movemap( core::kinematics::MoveMapCOP movemap) {

	using namespace core::kinematics;

	utility::vector1<core::Size> pivot_residues;

	for ( auto it=movemap->movemap_torsion_id_begin(), it_end=movemap->movemap_torsion_id_end(); it !=it_end; ++it ) {
		//Scaffold to new MM
		if ( it->first.second == core::id::BB && it->second ) {
			pivot_residues.push_back(it->first.first);
		}
	}
	return pivot_residues;
}

} //backrub
} //protocols

