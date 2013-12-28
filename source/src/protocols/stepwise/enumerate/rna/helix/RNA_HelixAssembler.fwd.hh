// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//  vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/Pose.fwd.hh
/// @brief  Pose forward declarations header
/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_HelixAssembler_fwd_hh
#define INCLUDED_protocols_rna_RNA_HelixAssembler_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace rna {
namespace helix {

	class RNA_HelixAssembler;
	typedef utility::pointer::owning_ptr< RNA_HelixAssembler > RNA_HelixAssemblerOP;

} //helix
} //rna
} //enumerate
} //stepwise
} //protocols

#endif
