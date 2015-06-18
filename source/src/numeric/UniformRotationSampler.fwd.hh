// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/UniformRotationSampler.fwd.hh
/// @brief  numeric::UniformRotationSampler forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_UniformRotationSampler_fwd_hh
#define INCLUDED_numeric_UniformRotationSampler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace numeric {


// Forward
class UniformRotationSampler;
typedef utility::pointer::shared_ptr< UniformRotationSampler > UniformRotationSamplerOP;
typedef utility::pointer::shared_ptr< UniformRotationSampler const > UniformRotationSamplerCOP;

} // namespace numeric


#endif // INCLUDED_numeric_UniformRotationSampler_FWD_HH

