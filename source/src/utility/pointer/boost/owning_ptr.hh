// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/boost/owning_ptr.hh
/// @brief  Owning smart pointer using boost header-only libraries
/// @author Luki Goldschmidt <lugo@uw.edu>


#ifndef INCLUDED_utility_pointer_boost_owning_ptr_hh
#define INCLUDED_utility_pointer_boost_owning_ptr_hh

#include <utility/pointer/boost/owning_ptr.fwd.hh>

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <utility/pointer/ReferenceCount.hh>

namespace utility {
namespace pointer {

	using boost::dynamic_pointer_cast;
	using boost::static_pointer_cast;
	using boost::const_pointer_cast;
	using boost::enable_shared_from_this;

}
}

#endif // INCLUDED_utility_pointer_boost_owning_ptr_hh
