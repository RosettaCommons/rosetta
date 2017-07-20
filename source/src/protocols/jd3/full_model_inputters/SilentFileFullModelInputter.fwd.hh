// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/SilentFileFullModelInputter.fwd.hh
/// @brief
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_full_model_inputters_SilentFileFullModelInputter_fwd_hh
#define INCLUDED_protocols_jd3_full_model_inputters_SilentFileFullModelInputter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

class SilentFileFullModelInputter;
typedef utility::pointer::shared_ptr< SilentFileFullModelInputter > SilentFileFullModelInputterOP;
typedef utility::pointer::shared_ptr< SilentFileFullModelInputter const> SilentFileFullModelInputterCOP;

} // full_model_inputters
} // jd3
} // protocols

#endif //INCLUDED_protocols_jd3_full_model_inputters_SilentFileFullModelInputter_FWD_HH
