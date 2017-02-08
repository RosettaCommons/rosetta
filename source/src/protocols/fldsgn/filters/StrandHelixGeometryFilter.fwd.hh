// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/filters/StrandHelixGeometryFilter.fwd.hh
/// @brief Another filter used in Marcos & Basanta et al. 2017 that needs to be updated.
/// @author Benjamin Basanta (basantab@uw.edu)

#ifndef INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_fwd_hh
#define INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace fldsgn {
namespace filters {

class StrandHelixGeometryFilter;

typedef utility::pointer::shared_ptr< StrandHelixGeometryFilter > StrandHelixGeometryFilterOP;
typedef utility::pointer::shared_ptr< StrandHelixGeometryFilter const > StrandHelixGeometryFilterCOP;

} //protocols
} //fldsgn
} //filters

#endif //INCLUDED_protocols_fldsgn_filters_StrandHelixGeometryFilter_fwd_hh
