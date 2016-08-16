// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EnergyMethod.fwd.hh
/// @brief  Energy Method class forward declaration and EnergyMethodType enum definition
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_methods_EnergyMethod_fwd_hh
#define INCLUDED_core_scoring_methods_EnergyMethod_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

enum EnergyMethodType {
	ci_2b = 1,
	cd_2b,
	ci_lr_2b,
	cd_lr_2b,
	ci_1b,
	cd_1b,
	ws,
	n_energy_method_types = ws // keep this guy last and equal to last type!
};

/// base class for the energy method hierarchy
class EnergyMethod;

typedef utility::pointer::shared_ptr< EnergyMethod > EnergyMethodOP;
typedef utility::pointer::shared_ptr< EnergyMethod const > EnergyMethodCOP;

}
}
}

#endif
