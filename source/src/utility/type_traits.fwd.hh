// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/type_traits.fwd.hh
/// @brief  Implemention of type traits related to Rosetta types, forward declarations
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_type_traits_fwd_hh
#define INCLUDED_utility_type_traits_fwd_hh

namespace utility {

#ifdef CXX11
template<typename T>
struct has_insertion_operator_s;
#endif

} // namespace utility


#endif // INCLUDED_utility_type_traits_fwd_hh
