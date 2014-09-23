// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/all.fwd.hh
/// @brief  utility::pointer package forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_pointer_all_fwd_hh
#define INCLUDED_utility_pointer_all_fwd_hh


// Classes
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

#ifdef PTR_REFCOUNT
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCountMI.fwd.hh>
#include <utility/pointer/ReferenceCountMI_.fwd.hh>
#endif

#endif // INCLUDED_utility_pointer_all_FWD_HH
