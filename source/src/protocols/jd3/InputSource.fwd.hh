// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/InputSource.fwd.hh
/// @brief  Forward declaration of the %InputSource class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_protocols_jd3_InputSource_HH
#define INCLUDED_protocols_jd3_InputSource_HH

//unit headers
#include <protocols/jd3/InputSource.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace jd3 {

class InputSource;

typedef utility::pointer::shared_ptr< InputSource > InputSourceOP;
typedef utility::pointer::shared_ptr< InputSource const > InputSourceCOP;

typedef utility::vector1< InputSourceOP > InputSources;

} // namespace jd3
} // namespace protocols

#endif //INCLUDED_protocols_jd3_InputSource_HH
