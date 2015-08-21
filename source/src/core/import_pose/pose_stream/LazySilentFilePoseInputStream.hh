// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/import_pose/pose_stream/LazySilentFilePoseInputStream.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_pose_stream_LazySilentFilePoseInputStream_HH
#define INCLUDED_core_import_pose_pose_stream_LazySilentFilePoseInputStream_HH

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/import_pose/pose_stream/LazySilentFilePoseInputStream.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <utility/file/FileName.hh>

#include <string>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {
namespace pose_stream {

class LazySilentFilePoseInputStream : public PoseInputStream {
	typedef std::string string;
	typedef utility::file::FileName FileName;

public:
	LazySilentFilePoseInputStream( utility::vector1< FileName > fns ) {
		filenames_        = fns;
		current_filename_ = filenames_.begin();

		sfd_.read_file( *current_filename_ );
		current_struct_   = sfd_.begin();
	}

	LazySilentFilePoseInputStream() {
		reset();
	}
	~LazySilentFilePoseInputStream() {}

	utility::vector1< FileName > const & filenames() const;

public:
	virtual bool has_another_pose();

	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set
	);
	virtual void fill_pose( core::pose::Pose& );

	virtual void reset();

private:
	utility::vector1< FileName > filenames_;
	utility::vector1< FileName >::const_iterator current_filename_;

	core::io::silent::SilentFileData sfd_;
	core::io::silent::SilentFileData::iterator current_struct_;
}; // LazySilentFilePoseInputStream

} // pose_stream
} // import_pose
} // core

#endif
