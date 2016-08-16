// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief   Forward declarations for MatchScoreWriter.
/// @author  Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_match_output_MatchScoreWriter_FWD_HH
#define INCLUDED_protocols_match_output_MatchScoreWriter_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace output {

/// @brief Simple class to write out a score file for the matcher
class MatchScoreWriter;

typedef utility::pointer::shared_ptr<MatchScoreWriter> MatchScoreWriterOP;
typedef utility::pointer::shared_ptr<MatchScoreWriter const> MatchScoreWriterCOP;

}  // namespace output
}  // namespace match
}  // namespace protocols

#endif  // INCLUDED_protocols_match_output_MatchScoreWriter_FWD_HH
