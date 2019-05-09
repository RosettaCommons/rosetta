// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/ImplicitLipidInfo.fwd.h
/// @brief An Implicit Lipid Membrane Model
/// @version menv-franklin2018
///
/// @detail This class defines physical and chemical properties of the membrane environment. It
/// is mainly used by the membrane energy function.
///  1. Parameters of the hydration function (e.g. thickness, rate, of transition, pore size)
///  2. Lipid composition details
///  3. Hydration function smoothing parameters
///  4. Structure-based lipid accessibiity information
///
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_ImplicitLipidInfo_fwd_hh
#define INCLUDED_core_conformation_membrane_ImplicitLipidInfo_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace core {
namespace conformation {
namespace membrane {

class ImplicitLipidInfo;

typedef utility::pointer::shared_ptr< ImplicitLipidInfo > ImplicitLipidInfoOP;
typedef utility::pointer::shared_ptr< ImplicitLipidInfo const > ImplicitLipidInfoCOP;

} //core
} //conformation
} //membrane

#endif //INCLUDED_core_conformation_membrane_ImplicitLipidInfo_fwd_hh
