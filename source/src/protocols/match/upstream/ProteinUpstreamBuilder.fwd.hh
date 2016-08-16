// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/upstream/ProteinUpstreamBuilder.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_ProteinUpstreamBuilder_fwd_hh
#define INCLUDED_protocols_match_upstream_ProteinUpstreamBuilder_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace upstream {

enum ChiStrategy {
	no_samples,   /// In particular, for proton chi, leave them fixed at some arbitrary value
	follow_EX_flags,
	rotameric_chi_mimic_EX_flags,
	rotameric_chi_step_by_value,
	rotameric_chi_step_wi_sd_range,
	rotameric_chi_partition_sd_range,
	nonrotameric_chi_sample_wi_nrchi_bin,
	nonrotameric_chi_sample_wi_nrchi_bin_to_lower_boundary
};

class SampleStrategyData;

class UpstreamResTypeGeometry;

typedef utility::pointer::shared_ptr< UpstreamResTypeGeometry > UpstreamResTypeGeometryOP;
typedef utility::pointer::shared_ptr< UpstreamResTypeGeometry const > UpstreamResTypeGeometryCOP;

class BuildSet;


class ProteinUpstreamBuilder;
typedef utility::pointer::shared_ptr< ProteinUpstreamBuilder > ProteinUpstreamBuilderOP;
typedef utility::pointer::shared_ptr< ProteinUpstreamBuilder const > ProteinUpstreamBuilderCOP;


}
}
}

#endif

