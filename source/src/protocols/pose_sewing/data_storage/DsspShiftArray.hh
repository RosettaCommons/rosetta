// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/DsspShiftArray.hh
/// @brief A container for an array of secondary structure boundaries
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_pose_sewing_data_storage_DsspShiftArray_hh
#define INCLUDED_protocols_pose_sewing_data_storage_DsspShiftArray_hh

#include <protocols/pose_sewing/data_storage/DsspShiftArray.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>

namespace protocols {
namespace pose_sewing {
namespace data_storage {

/// @brief A container for an array of secondary structure boundaries
class DsspShiftArray : public utility::VirtualBase {

public:

	DsspShiftArray();

	DsspShiftArray(core::pose::Pose const & pose);

	DsspShiftArray(DsspShiftArray const & src);

	virtual ~DsspShiftArray();

	DsspShiftArrayOP
	clone() const;

	void
	populate(core::pose::Pose const & pose);

	core::Size
	get_distance_to_nth_shift(core::Size N_res, core::Size shift_gap);

private:
	utility::vector1<core::Size> shift_array_;
};


} //protocols
} //pose_sewing
} //data_storage



#endif //INCLUDED_protocols_pose_sewing_data_storage_DsspShiftArray_hh





