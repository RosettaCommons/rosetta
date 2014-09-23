// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/ProcessorFactory.hh
/// @brief  Declaration for a factory class responsible for instantiating
///         MatchProcessor classes.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_ProcessorFactory_fwd_hh
#define INCLUDED_protocols_match_output_ProcessorFactory_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace output {

class ProcessorFactory;
typedef utility::pointer::shared_ptr< ProcessorFactory > ProcessorFactoryOP;
typedef utility::pointer::shared_ptr< ProcessorFactory const > ProcessorFactoryCOP;

}
}
}

#endif
