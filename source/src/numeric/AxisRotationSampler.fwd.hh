// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/AxisRotationSampler.fwd.hh
/// @brief  numeric::AxisRotationSampler forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_AxisRotationSampler_fwd_hh
#define INCLUDED_numeric_AxisRotationSampler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace numeric {


// Forward
class AxisRotationSampler;
typedef utility::pointer::shared_ptr< AxisRotationSampler > AxisRotationSamplerOP;
typedef utility::pointer::shared_ptr< AxisRotationSampler const > AxisRotationSamplerCOP;

} // namespace numeric


#endif // INCLUDED_numeric_AxisRotationSampler_FWD_HH

