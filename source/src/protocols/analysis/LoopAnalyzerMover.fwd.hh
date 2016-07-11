// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/analysis/LoopAnalyzerMover.fwd.hh
/// @brief  fwd declaration for LoopAnalyzer: intense analysis of loop quality
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_analysis_LoopAnalyzerMover_fwd_HH
#define INCLUDED_protocols_analysis_LoopAnalyzerMover_fwd_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace analysis {

//Forwards and OP typedefs
class LoopAnalyzerMover;
typedef utility::pointer::shared_ptr< LoopAnalyzerMover > LoopAnalyzerMoverOP;

}//analysis
}//protocols

#endif //INCLUDED_protocols_analysis_LoopAnalyzerMover_FWD_HH
