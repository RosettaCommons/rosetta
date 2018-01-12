// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/FullModelInputSource.fwd.hh
/// @brief  Forward declaration of the %FullModelInputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_full_model_inputters_FullModelInputSource_HH
#define INCLUDED_protocols_jd3_full_model_inputters_FullModelInputSource_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace jd3 {
namespace full_model_inputters {

class FullModelInputSource;

typedef utility::pointer::shared_ptr< FullModelInputSource > FullModelInputSourceOP;
typedef utility::pointer::shared_ptr< FullModelInputSource const > FullModelInputSourceCOP;

typedef utility::vector1< FullModelInputSourceOP > FullModelInputSources;

} // namespace full_model_inputters
} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_FullModelInputSource_HH
