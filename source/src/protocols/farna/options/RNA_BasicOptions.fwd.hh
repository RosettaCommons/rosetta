// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/options/RNA_BasicOptions.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_BasicOptions_FWD_HH
#define INCLUDED_protocols_farna_RNA_BasicOptions_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace farna {
namespace options {
	
	class RNA_BasicOptions;
	typedef utility::pointer::shared_ptr< RNA_BasicOptions > RNA_BasicOptionsOP;
	typedef utility::pointer::shared_ptr< RNA_BasicOptions const > RNA_BasicOptionsCOP;
	
} //options
} //farna
} //protocols

#endif
