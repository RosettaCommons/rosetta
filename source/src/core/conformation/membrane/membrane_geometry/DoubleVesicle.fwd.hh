// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/membrane_geometry/DoubleVesicle.fwd.hh
/// @brief    Data describing the parameters of a vesicle
///
/// @details  DoubleVesicle class contains the parameters of a vesicle
///      and the function to calculate the transition from the
///      hydrophobic environment of the bilayer to a hydrophilic environment.
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_geometry().
///
///
/// @author   Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_fwd_hh
#define INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {
namespace membrane_geometry {

/// @brief DoubleVesicle describes a flat membrane extending into the xy plane
class DoubleVesicle;
typedef utility::pointer::shared_ptr< DoubleVesicle > DoubleVesicleOP;
typedef utility::pointer::shared_ptr< DoubleVesicle const > DoubleVesicleCOP;

} // membrane_geometry
} // membrane
} // conformation
} // core

#endif //INCLUDED_core_conformation_membrane_geometry_DoubleVesicle_fwd_hh

