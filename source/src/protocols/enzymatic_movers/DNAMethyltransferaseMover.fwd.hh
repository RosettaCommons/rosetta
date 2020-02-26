// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzymatic_movers/DNAMethyltransferaseMover.fwd.hh
/// @brief  Forward declarations for DNAMethyltransferaseMover
/// @author Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_protocols_enzymatic_movers_DNAMethyltransferaseMover_FWD_HH
#define INCLUDED_protocols_enzymatic_movers_DNAMethyltransferaseMover_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace enzymatic_movers {

/// @brief  An EnzymaticMover class for simulating the DNA-methylation of a pose.
class DNAMethyltransferaseMover;

typedef utility::pointer::shared_ptr< DNAMethyltransferaseMover > DNAMethyltransferaseMoverOP;
typedef utility::pointer::shared_ptr< DNAMethyltransferaseMover const > DNAMethyltransferaseMoverCOP;

}  // namespace enzymatic_movers
}  // namespace protocols

#endif  // INCLUDED_protocols_enzymatic_movers_DNAMethyltransferaseMover_FWD_HH
