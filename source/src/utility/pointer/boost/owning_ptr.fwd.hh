// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/boost/owning_ptr.hh
/// @brief  Owning smart pointer using C++11
/// @author Luki Goldschmidt <lugo@uw.edu>


#ifndef INCLUDED_utility_pointer_boost_owning_ptr_fwd_hh
#define INCLUDED_utility_pointer_boost_owning_ptr_fwd_hh

#include <boost/shared_ptr.hpp>

namespace utility {
namespace pointer {

using boost::shared_ptr;

}
}

#endif // INCLUDED_utility_pointer_boost_owning_ptr_fwd_hh
