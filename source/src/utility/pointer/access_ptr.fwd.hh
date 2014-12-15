// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/access_ptr.fwd.hh
/// @brief  Non-owning access smart pointer -- dispatch class
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifdef PTR_STD
#include <utility/pointer/std/access_ptr.hh>
#endif

#ifdef PTR_BOOST
#include <utility/pointer/boost/access_ptr.hh>
#endif
