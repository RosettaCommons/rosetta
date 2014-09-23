// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/ExternalGeomSampler.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_ExternalGeomSampler_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_ExternalGeomSampler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

enum ExternalTransform {
	HT_tor_U3D1 = 1,
	HT_ang_U2D1,
	HT_tor_U2D2,
	HT_ang_U1D2,
	HT_tor_U1D3,
	n_external_transforms = HT_tor_U1D3
};

class ExternalGeomSampler;

typedef utility::pointer::shared_ptr< ExternalGeomSampler > ExternalGeomSamplerOP;
typedef utility::pointer::shared_ptr< ExternalGeomSampler const > ExternalGeomSamplerCOP;

}
}
}

#endif
