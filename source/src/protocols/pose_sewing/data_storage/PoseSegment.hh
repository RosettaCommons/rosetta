// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/PoseSegment.hh
/// @brief a region of a Pose with contiguous secondary structure
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_pose_sewing_data_storage_PoseSegment_hh
#define INCLUDED_protocols_pose_sewing_data_storage_PoseSegment_hh

#include <protocols/pose_sewing/data_storage/PoseSegment.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace pose_sewing {
namespace data_storage {

/// @brief a region of a Pose with contiguous secondary structure
class PoseSegment : public utility::VirtualBase {

public:

	PoseSegment();
	PoseSegment(core::Size,core::Size,char,core::pose::PoseCOP);
	PoseSegment(PoseSegment const & src);
	PoseSegment( PoseSegment const & src, core::Size start,core::Size end);

	~PoseSegment() override;

	PoseSegmentOP
	clone() const;
	void
	set_starting_residue(core::Size new_starting);
	void
	set_ending_residue(core::Size new_ending);
	void
	set_dssp_code(char new_code);
	void
	set_source_pose(core::pose::PoseCOP new_posecop);

	core::Size
	get_starting_residue() const;
	core::Size
	get_ending_residue() const;
	char
	get_dssp_code() const;
	core::pose::PoseCOP
	get_source_pose() const;
	core::Size
	get_length() const;
private:
	core::Size starting_residue_;
	core::Size ending_residue_;
	char dssp_code_;
	core::pose::PoseCOP source_pose_;
};


} //protocols
} //pose_sewing
} //data_storage



#endif //INCLUDED_protocols_pose_sewing_data_storage_PoseSegment_hh





