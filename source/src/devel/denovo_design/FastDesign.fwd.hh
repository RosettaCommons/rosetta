// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/denovo_design/FastDesign.fwd.hh
/// @brief  FastDesign forward header
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_FastDesign_fwd_hh
#define INCLUDED_devel_denovo_design_FastDesign_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace devel {
namespace denovo_design {

// Forward
class FastDesign;

// Types
typedef  utility::pointer::shared_ptr< FastDesign >  FastDesignOP;
typedef  utility::pointer::shared_ptr< FastDesign const >  FastDesignCOP;

typedef  utility::pointer::weak_ptr< FastDesign >  FastDesignAP;
typedef  utility::pointer::weak_ptr< FastDesign const >  FastDesignCAP;


} // namespace denovo_design
} // namespace devel

#endif
