// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/ProteinSCSampler.fwd.hh
/// @brief  Forward declaration for base class and simple derived class for SC sampling in matcher
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_ProteinSCSampler_fwd_hh
#define INCLUDED_protocols_match_upstream_ProteinSCSampler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace upstream {

class ProteinSCSampler;

typedef utility::pointer::shared_ptr< ProteinSCSampler > ProteinSCSamplerOP;
typedef utility::pointer::shared_ptr< ProteinSCSampler const > ProteinSCSamplerCOP;

class DunbrackSCSampler;

typedef utility::pointer::shared_ptr< DunbrackSCSampler > DunbrackSCSamplerOP;
typedef utility::pointer::shared_ptr< DunbrackSCSampler const > DunbrackSCSamplerCOP;

}
}
}

#endif
