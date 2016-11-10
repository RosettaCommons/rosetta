// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/MPLipidAccessibility.fwd.hh
/// @brief Mover that computes which residues are lipid accessible and puts that information into the B-factors: 50 is lipid accessible, 0 is lipid INaccessible
/// @author Julia Koehler Leman (julia.koehler1982@gmail.com)


#ifndef INCLUDED_protocols_membrane_MPLipidAccessibility_fwd_hh
#define INCLUDED_protocols_membrane_MPLipidAccessibility_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace membrane {

class MPLipidAccessibility;

typedef utility::pointer::shared_ptr< MPLipidAccessibility > MPLipidAccessibilityOP;
typedef utility::pointer::shared_ptr< MPLipidAccessibility const > MPLipidAccessibilityCOP;



} //protocols
} //membrane


#endif //INCLUDED_protocols_membrane_MPLipidAccessibility_fwd_hh





