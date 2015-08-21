// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/scoring/methods/saxs/PDDFEnergy.fwd.hh
/// @brief Forward declaration of a PDDF energy term
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_scoring_methods_saxs_PDDFEnergy_fwd_hh
#define INCLUDED_protocols_scoring_methods_saxs_PDDFEnergy_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace scoring {
namespace methods {
namespace saxs {

class PDDFEnergy;

typedef utility::pointer::shared_ptr<PDDFEnergy>
	PDDFEnergyOP;
typedef utility::pointer::shared_ptr<PDDFEnergy const>
	PDDFEnergyCOP;

} // saxs
} // methods
} // scoring
} //protocols

#endif /* INCLUDED_protocols_scoring_methods_saxs_PDDFEnergy_FWD_HH */
