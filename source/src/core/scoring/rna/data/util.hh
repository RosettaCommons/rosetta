// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/rna/data/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_data_util_HH
#define INCLUDED_core_scoring_rna_data_util_HH

#include <core/types.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace scoring {
namespace rna {
namespace data {

core::Size
lookup_idx( core::Real const value, utility::vector1< core::Real > & values );

core::Size
get_bool_idx( bool const value, utility::vector1< bool > const & values );

core::Size
get_idx( core::Real const value, utility::vector1< core::Real > const & values );

} //data
} //rna
} //scoring
} //core

#endif
