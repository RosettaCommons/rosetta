// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/denovo_design/filters/SSShapeComplementarityFilter.fwd.hh
/// @brief  Fwd declarations for Tom's denovo design protocol
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_filters_SSShapeComplementarityFilter_fwd_hh
#define INCLUDED_protocols_denovo_design_filters_SSShapeComplementarityFilter_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {
namespace filters {
// Forward
class SSShapeComplementarityFilter;

// Types
typedef  utility::pointer::owning_ptr< SSShapeComplementarityFilter >  SSShapeComplementarityFilterOP;
typedef  utility::pointer::owning_ptr< SSShapeComplementarityFilter const >  SSShapeComplementarityFilterCOP;

typedef  utility::pointer::access_ptr< SSShapeComplementarityFilter >  SSShapeComplementarityFilterAP;
typedef  utility::pointer::access_ptr< SSShapeComplementarityFilter const >  SSShapeComplementarityFilterCAP;

} // namespace filters
} // namespace kinematics
} // namespace core

#endif
