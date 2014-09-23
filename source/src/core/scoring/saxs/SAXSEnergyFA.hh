// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/SAXSEnergy.hh
/// @brief  "Energy" based on a similarity of theoretical SAXS spectrum computed for a pose and the experimental data
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_core_scoring_saxs_SAXSEnergyFA_hh
#define INCLUDED_core_scoring_saxs_SAXSEnergyFA_hh

// Package headers
// AUTO-REMOVED #include <core/scoring/saxs/FormFactorManager.hh>
#include <core/scoring/saxs/FormFactor.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>
#include <core/scoring/saxs/SAXSEnergy.hh>

// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <string>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.fwd.hh>
// AUTO-REMOVED #include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
// AUTO-REMOVED #include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/SecondaryStructureWeights.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
// AUTO-REMOVED #include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/scoring/saxs/DistanceHistogram.fwd.hh>
#include <core/scoring/saxs/DistanceHistogram.hh>
#include <core/scoring/saxs/FormFactor.fwd.hh>
#include <core/scoring/saxs/FormFactorManager.fwd.hh>
#include <core/scoring/saxs/SAXSEnergyCreator.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorCEN.hh>
#include <utility/down_cast.hh>
// AUTO-REMOVED #include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <numeric/sphericalVector.fwd.hh>
// AUTO-REMOVED #include <numeric/trig.functions.hh>
#include <numeric/types.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>
#include <numeric/interpolation/spline/Interpolator.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <algorithm>
#include <cassert>
// AUTO-REMOVED #include <cmath>
#include <cstddef>
// AUTO-REMOVED #include <cstdlib>
// AUTO-REMOVED #include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
// AUTO-REMOVED #include <set>
// AUTO-REMOVED #include <sstream>
#include <vector>


namespace core {
namespace scoring {
namespace saxs {


class SAXSEnergyFA : public SAXSEnergy  {
public:

    /// c-tor
    SAXSEnergyFA() : SAXSEnergy( fa_cfg_file_,
	    chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD),saxs_fa_score,methods::EnergyMethodCreatorOP( new SAXSEnergyCreatorFA )) {}

};


}
}
}

#endif
