// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/SegmentEnvelope.fwd.hh
/// @brief a descriptor of properties a SewAnything secondar structral segment may have
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_fwd_hh
#define INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace pose_sewing {
namespace data_storage {

class SegmentEnvelope;

typedef utility::pointer::shared_ptr< SegmentEnvelope > SegmentEnvelopeOP;
typedef utility::pointer::shared_ptr< SegmentEnvelope const > SegmentEnvelopeCOP;

} //protocols
} //pose_sewing
} //data_storage

#endif //INCLUDED_protocols_pose_sewing_data_storage_SegmentEnvelope_fwd_hh
