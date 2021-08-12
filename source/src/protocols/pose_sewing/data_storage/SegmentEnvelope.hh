// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/SegmentEnvelope.hh
/// @brief a descriptor of properties a SewAnything secondar structral segment may have
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_hh
#define INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_hh

#include <protocols/pose_sewing/data_storage/SegmentEnvelope.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <core/types.hh>
#include <protocols/pose_sewing/data_storage/PoseSegment.hh>

namespace protocols {
namespace pose_sewing {
namespace data_storage {

/// @brief a descriptor of properties a SewAnything secondar structral segment may have
class SegmentEnvelope : public utility::VirtualBase {

public:

	SegmentEnvelope();
	SegmentEnvelope(SegmentEnvelope const & src);

	~SegmentEnvelope() override;

	SegmentEnvelopeOP
	clone() const;

	void
	set_minimum_length(core::Size);

	void
	set_maximum_length(core::Size);

	void
	set_permissible_secondary_structures(std::string &);

	core::Size
	get_minimum_length();

	core::Size
	get_maximum_length();

	std::string
	get_permissible_secondary_structures();

	bool
	is_permissible_secondary_structure(char);

	bool
	is_permissible_length(core::Size);

	bool
	is_valid(char,core::Size);

	bool
	is_valid(PoseSegmentOP);

private:

	core::Size minimum_length_;
	core::Size maximum_length_;
	std::string permissible_secondary_structures_;
};


} //protocols
} //pose_sewing
} //data_storage



#endif //INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_hh





