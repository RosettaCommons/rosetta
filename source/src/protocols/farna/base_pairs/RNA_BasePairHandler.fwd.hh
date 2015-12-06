// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/base_pairs/RNA_BasePairHandler.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_BasePairHandler_FWD_HH
#define INCLUDED_protocols_farna_RNA_BasePairHandler_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace farna {
namespace base_pairs {
	
	class RNA_BasePairHandler;
	typedef utility::pointer::shared_ptr< RNA_BasePairHandler > RNA_BasePairHandlerOP;
	typedef utility::pointer::shared_ptr< RNA_BasePairHandler const > RNA_BasePairHandlerCOP;
	
} //base_pairs
} //farna
} //protocols

#endif
