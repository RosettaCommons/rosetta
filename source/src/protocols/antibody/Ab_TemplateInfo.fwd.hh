// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocls/antibody/Ab_TemplateInfo.fwd.hh
/// @brief grafts a cdr onto the template of an antibody framework
/// @details
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_Ab_TemplateInfo_fwd_hh
#define INCLUDED_protocols_antibody_Ab_TemplateInfo_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {

// Forward
class Ab_TemplateInfo;

typedef utility::pointer::shared_ptr< Ab_TemplateInfo > Ab_TemplateInfoOP;
typedef utility::pointer::shared_ptr< Ab_TemplateInfo const > Ab_TemplateInfoCOP;

} // antibody
} // protocols

#endif


