// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/denovo_design/BuildSheet.fwd.hh
/// @brief  BuildSheet forward header
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_BuildSheet_fwd_hh
#define INCLUDED_devel_denovo_design_BuildSheet_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace devel {
namespace denovo_design {

// Forward
class BuildSheet;

// Types
typedef  utility::pointer::shared_ptr< BuildSheet >  BuildSheetOP;
typedef  utility::pointer::shared_ptr< BuildSheet const >  BuildSheetCOP;

typedef  utility::pointer::weak_ptr< BuildSheet >  BuildSheetAP;
typedef  utility::pointer::weak_ptr< BuildSheet const >  BuildSheetCAP;


} // namespace denovo_design
} // namespace devel

#endif
