// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/import_pose/pose_stream/PDBPoseInputStream.cc
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>

#include <basic/datacache/BasicDataCache.hh>

#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

#include <basic/datacache/CacheableString.hh>

// C++ headers
#include <string>
#include <utility/exit.hh>

#include <utility/vector1.hh>

namespace core {
namespace import_pose {
namespace pose_stream {


void PDBPoseInputStream::set_filenames(
	utility::vector1< utility::file::FileName > filenames
) {
	filenames_ = filenames;
	current_position_ = filenames_.begin();
}

utility::vector1< utility::file::FileName > PDBPoseInputStream::get_filenames() {
	return filenames_;
}

bool PDBPoseInputStream::has_another_pose() {
	return ( current_position_ != filenames_.end() );
}

void PDBPoseInputStream::reset(){
	current_position_ = filenames_.begin();
}

void PDBPoseInputStream::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set,
	bool const //metapatches /*= true*/
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"PDBPoseInputStream: called fill_pose, but I have no more Poses!"
		);
	}

	core::import_pose::pose_from_file( pose, residue_set, *current_position_ , core::import_pose::PDB_file);

	// set up a tag using input filename.
	using namespace basic::datacache;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( *current_position_ ) )
	);
	++current_position_;
	preprocess_pose( pose );
}

void PDBPoseInputStream::fill_pose(
	core::pose::Pose & pose,
	bool const //metapatches /*= true*/
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"PDBPoseInputStream: called fill_pose, but I have no more Poses!"
		);
	}

	core::import_pose::pose_from_file( pose, *current_position_ , core::import_pose::PDB_file);

	// set up a tag using input filename.
	using namespace basic::datacache;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( *current_position_ ) )
	);
	++current_position_;
	preprocess_pose( pose );
}

utility::vector1< core::pose::PoseOP > PDBPoseInputStream::get_all_poses(
	core::chemical::ResidueTypeSet const & residue_set
) {
	utility::vector1< core::pose::PoseOP > pose_list;
	pose_list.resize( filenames_.size() );

	while ( has_another_pose() ) {
		core::pose::PoseOP pose( new core::pose::Pose );
		fill_pose( *pose, residue_set );
		pose_list.push_back( pose );
	}

	return pose_list;
}

/// @brief adds a list of files each containing lists of PDBs
void PDBPoseInputStream::add_list_filenames(
	utility::vector1< utility::file::FileName > list_fns
) {
	using utility::vector1;
	bool init_current_position( filenames_.size() == 0 );

	for ( utility::file::FileName const & filename_obj : list_fns ) {
		std::string const & filename( filename_obj.name() );
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}

		std::string line;
		while ( getline(data, line) ) {
			filenames_.push_back( utility::file::FileName(line) );
		}
		data.close();
	}

	if ( init_current_position ) {
		current_position_ = filenames_.begin();
	}
} // add_list_files

} // pose_stream
} // import_pose
} // core
