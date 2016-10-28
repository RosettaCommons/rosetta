// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/GlycosyltransferaseMover.fwd.hh
/// @brief  Forward declarations for GlycosyltransferaseMover
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_FWD_HH
#define INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace enzymatic_movers {

/// @brief  An EnzymaticMover class for simulating the glycosylation of a pose.
class GlycosyltransferaseMover;

typedef utility::pointer::shared_ptr< GlycosyltransferaseMover > GlycosyltransferaseMoverOP;
typedef utility::pointer::shared_ptr< GlycosyltransferaseMover const > GlycosyltransferaseMoverCOP;

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_GlycosyltransferaseMover_FWD_HH
