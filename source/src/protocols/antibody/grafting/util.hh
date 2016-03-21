// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/util.hh
/// @brief Helpers for antibody grafting code
/// @author Sergey Lyskov

#ifndef INCLUDED_protocols_antibody_grafting_util_hh
#define INCLUDED_protocols_antibody_grafting_util_hh

#ifdef CXX11

#ifdef __clang__

		#define __ANTIBODY_GRAFTING__
#else
#if (__GNUC__ > 3  &&  __GNUC_MINOR__ > 8)  || (__GNUC__ > 4) // We need at least GCC-4.9 to compiler Antibody code
#define __ANTIBODY_GRAFTING__
#endif
#endif

#endif // CXX11


namespace protocols {
namespace antibody {
namespace grafting {

/// Check if regex library is functional
bool antibody_grafting_usable();


} // grafting
} // antibody
} // protocols

#endif // __ANTIBODY_GRAFTING__




#ifdef __ANTIBODY_GRAFTING__

namespace protocols {
namespace antibody {
namespace grafting {

} // grafting
} // antibody
} // protocols

#endif // INCLUDED_protocols_antibody_grafting_util_hh
