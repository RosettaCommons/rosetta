// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/AqueousPoreParameters.fwd.hh
/// @brief A class for defining an aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_AqueousPoreParameters_fwd_hh
#define INCLUDED_core_conformation_membrane_AqueousPoreParameters_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace conformation {
namespace membrane {

class AqueousPoreParameters;

typedef utility::pointer::shared_ptr< AqueousPoreParameters > AqueousPoreParametersOP;
typedef utility::pointer::shared_ptr< AqueousPoreParameters const > AqueousPoreParametersCOP;

} //core
} //conformation
} //membrane

#endif //INCLUDED_core_conformation_membrane_AqueousPoreParameters_fwd_hh
