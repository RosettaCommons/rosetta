// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/SS_Info.fwd.hh
/// @brief  Data cache forward declarations for the secondary-structure based scores.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_SS_Info_fwd_hh
#define INCLUDED_core_scoring_SS_Info_fwd_hh

/// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class BB_Pos;
struct Strands;
struct Helices;
class SS_Info;

typedef utility::pointer::shared_ptr< SS_Info > SS_InfoOP;
typedef utility::pointer::shared_ptr< SS_Info const > SS_InfoCOP;

} // namespace scoring
} // namespace core

#endif
