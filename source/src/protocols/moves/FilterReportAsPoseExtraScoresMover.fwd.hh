// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/FilterReportAsPoseExtraScoresMover.fwd.hh
/// @brief This Mover runs a Filter and dumps the report_sm value to Pose's extra scores (for later JD reporting)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_moves_FilterReportAsPoseExtraScoresMover_fwd_hh
#define INCLUDED_protocols_moves_FilterReportAsPoseExtraScoresMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace moves {

class FilterReportAsPoseExtraScoresMover;

typedef utility::pointer::shared_ptr< FilterReportAsPoseExtraScoresMover > FilterReportAsPoseExtraScoresMoverOP;
typedef utility::pointer::shared_ptr< FilterReportAsPoseExtraScoresMover const > FilterReportAsPoseExtraScoresMoverCOP;



} //protocols
} //moves


#endif //INCLUDED_protocols_moves_FilterReportAsPoseExtraScoresMover_fwd_hh





