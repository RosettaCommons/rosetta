// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/VectorPoseFilter.fwd.hh
/// @brief Used for filters that require multiple poses to act upon
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_PROTOCOLS_FILTERS_VectorPoseFilter_FWD_HH
#define INCLUDED_PROTOCOLS_FILTERS_VectorPoseFilter_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace filters {

class VectorPoseFilter;
typedef utility::pointer::shared_ptr<VectorPoseFilter> VectorPoseFilterOP;
typedef utility::pointer::shared_ptr<VectorPoseFilter const> VectorPoseFilterCOP;

}  // namespace filters
}  // namespace protocols

#endif  // INCLUDED_PROTOCOLS_FILTERS_VectorPoseFilter_FWD_HH
