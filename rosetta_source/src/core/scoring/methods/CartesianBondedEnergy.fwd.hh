// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.fwd.hh
/// @brief
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_methods_CartesianBondedEnergy_fwd_hh
#define INCLUDED_core_scoring_methods_CartesianBondedEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class IdealParametersDatabase;
typedef  utility::pointer::access_ptr< IdealParametersDatabase > IdealParametersDatabaseOPAP;
typedef  utility::pointer::access_ptr< IdealParametersDatabase const > IdealParametersDatabaseOPCAP;
typedef  utility::pointer::owning_ptr< IdealParametersDatabase > IdealParametersDatabaseOP;
typedef  utility::pointer::owning_ptr< IdealParametersDatabase const > IdealParametersDatabaseCOP;

class CartesianBondedEnergy;
typedef  utility::pointer::access_ptr< CartesianBondedEnergy > CartesianBondedEnergyAP;
typedef  utility::pointer::access_ptr< CartesianBondedEnergy const > CartesianBondedEnergyCAP;
typedef  utility::pointer::owning_ptr< CartesianBondedEnergy > CartesianBondedEnergyOP;
typedef  utility::pointer::owning_ptr< CartesianBondedEnergy const > CartesianBondedEnergyCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_CartesianBondedEnergy_HH
