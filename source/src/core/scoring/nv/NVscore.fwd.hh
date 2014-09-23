// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/NV/NVscore.fwd.hh
/// @brief  Typedefs and forward declarations for class NVscore
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nv_NVscore_fwd_hh
#define INCLUDED_core_scoring_nv_NVscore_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace nv {

class NVscore;

typedef utility::pointer::shared_ptr< NVscore > NVscoreOP;
typedef utility::pointer::shared_ptr< NVscore const > NVscoreCOP;

} //NV
} //scoring
} //core

#endif /* NVSCORE_FWD_HH_ */
