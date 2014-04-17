// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/full_model_info/FullModelInfoSetupFromCommandLine.fwd.hh
/// @brief 
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_full_model_info_FullModelInfoSetupFromCommandLine_FWD_HH
#define INCLUDED_core_pose_full_model_info_FullModelInfoSetupFromCommandLine_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace full_model_info {
	
	class FullModelInfoSetupFromCommandLine;
	typedef utility::pointer::owning_ptr< FullModelInfoSetupFromCommandLine > FullModelInfoSetupFromCommandLineOP;
	typedef utility::pointer::owning_ptr< FullModelInfoSetupFromCommandLine const > FullModelInfoSetupFromCommandLineCOP;
	
} //full_model_info 
} //pose 
} //core 

#endif
