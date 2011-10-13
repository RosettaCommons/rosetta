// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/filters/ShapeComplementarityFilter.fwd.hh
/// @brief
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)


#ifndef INCLUDED_protocols_filters_ShapeComplementarityFilter_fwd_hh
#define INCLUDED_protocols_filters_ShapeComplementarityFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace filters {

// Forward
class ShapeComplementarityFilter;

// Types
typedef utility::pointer::owning_ptr< ShapeComplementarityFilter >  ShapeComplementarityFilterOP;
typedef utility::pointer::owning_ptr< ShapeComplementarityFilter const >  ShapeComplementarityFilterCOP;

} // namespace protocols
} // namespace filters

#endif
