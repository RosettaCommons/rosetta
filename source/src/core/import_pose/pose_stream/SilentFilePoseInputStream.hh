// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/import_pose/pose_stream/SilentFilePoseInputStream.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_pose_stream_SilentFilePoseInputStream_HH
#define INCLUDED_core_import_pose_pose_stream_SilentFilePoseInputStream_HH

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentFileData.hh> // for its ::iterator
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <utility/file/FileName.hh> // because vector of FileName + destructor in .hh

#include <string>

#include <utility/vector1.fwd.hh>


namespace core {
namespace import_pose {
namespace pose_stream {


class SilentFilePoseInputStream : public PoseInputStream {

	typedef std::string string;
	typedef utility::file::FileName FileName;

	// constructors
public:
	SilentFilePoseInputStream();
	SilentFilePoseInputStream( utility::vector1< FileName > const & fns );
	SilentFilePoseInputStream( std::string const & fn );
	SilentFilePoseInputStream(
		utility::vector1< FileName > const & fns,
		bool order_by_energy
	);

	SilentFilePoseInputStream(
		utility::vector1< FileName > const & fns,
		core::Real energy_cut
	);

	SilentFilePoseInputStream(
		utility::vector1< FileName > const & fns,
		utility::vector1< string > const & input_tags
	);

	SilentFilePoseInputStream(
		utility::vector1< FileName > const & fns,
		utility::vector1< string > const & input_tags,
		core::Real energy_cut
	);

	void
	set_silent_file_data( core::io::silent::SilentFileDataOP & sfd );

	/// @brief Access the SilentFileData owning pointer directly.
	///
	core::io::silent::SilentFileDataOP silent_file_data() { return sfd_; }

	void
	set_record_source( bool const & setting );

	~SilentFilePoseInputStream() {}

public: // methods specific to SilentFilePoseInputStream class
	void renumber_decoys( bool const setting );
	void tags( utility::vector1< string > const & tags );
	void filenames( utility::vector1< FileName > const & filenames );

	core::Real energy_cut() const;
	bool renumber_decoys() const;
	utility::vector1< string > tags() const;
	utility::vector1< FileName > filenames() const;

	void set_order_by_energy( bool const & setting );

public: // class-wide methods
	virtual bool has_another_pose();

	virtual void reset();

	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set,
		bool const metapatches = true
	);
	virtual void fill_pose( core::pose::Pose&, bool const metapatches = true );

	core::io::silent::SilentStructOP next_struct();

private:
	utility::vector1< string > tags_;
	utility::vector1< FileName > filenames_;
	core::io::silent::SilentFileDataOP sfd_;
	core::io::silent::SilentFileData::iterator current_position_;
	bool renumber_decoys_;
	core::Real energy_cut_;
	bool order_by_energy_;
	bool record_source_;

	void read_all_files_();
}; // SilentFilePoseInputStream

} // pose_stream
} // import_pose
} // core

#endif
