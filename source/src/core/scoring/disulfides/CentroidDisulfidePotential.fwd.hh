// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/disulfides/CentroidDisulfidePotential.fwd.hh
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   1/27/09


#ifndef INCLUDED_core_scoring_disulfides_CentroidDisulfidePotential_fwd_hh
#define INCLUDED_core_scoring_disulfides_CentroidDisulfidePotential_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace disulfides {

class CentroidDisulfidePotential;
typedef utility::pointer::shared_ptr< CentroidDisulfidePotential > CentroidDisulfidePotentialOP;
typedef utility::pointer::shared_ptr< CentroidDisulfidePotential const > CentroidDisulfidePotentialCOP;
typedef utility::pointer::weak_ptr< CentroidDisulfidePotential > CentroidDisulfidePotentialAP;
typedef utility::pointer::weak_ptr< CentroidDisulfidePotential const > CentroidDisulfidePotentialCAP;

class Cb_Distance_Func;
class Cen_Distance_Func;
class CaCbCb_Angle_Func;
class NCaCaC_Dihedral_Func;
class CaCbCbCa_Dihedral_Func;

typedef utility::pointer::shared_ptr< Cb_Distance_Func > Cb_Distance_FuncOP;
typedef utility::pointer::shared_ptr< Cen_Distance_Func > Cen_Distance_FuncOP;
typedef utility::pointer::shared_ptr< CaCbCb_Angle_Func > CaCbCb_Angle_FuncOP;
typedef utility::pointer::shared_ptr< NCaCaC_Dihedral_Func > NCaCaC_Dihedral_FuncOP;
typedef utility::pointer::shared_ptr< CaCbCbCa_Dihedral_Func > CaCbCbCa_Dihedral_FuncOP;

typedef utility::pointer::shared_ptr< Cb_Distance_Func const > Cb_Distance_FuncCOP;
typedef utility::pointer::shared_ptr< Cen_Distance_Func const > Cen_Distance_FuncCOP;
typedef utility::pointer::shared_ptr< CaCbCb_Angle_Func const > CaCbCb_Angle_FuncCOP;
typedef utility::pointer::shared_ptr< NCaCaC_Dihedral_Func const > NCaCaC_Dihedral_FuncCOP;
typedef utility::pointer::shared_ptr< CaCbCbCa_Dihedral_Func const > CaCbCbCa_Dihedral_FuncCOP;

typedef utility::pointer::weak_ptr< Cb_Distance_Func > Cb_Distance_FuncAP;
typedef utility::pointer::weak_ptr< Cen_Distance_Func > Cen_Distance_FuncAP;
typedef utility::pointer::weak_ptr< CaCbCb_Angle_Func > CaCbCb_Angle_FuncAP;
typedef utility::pointer::weak_ptr< NCaCaC_Dihedral_Func > NCaCaC_Dihedral_FuncAP;
typedef utility::pointer::weak_ptr< CaCbCbCa_Dihedral_Func > CaCbCbCa_Dihedral_FuncAP;

typedef utility::pointer::weak_ptr< Cb_Distance_Func const > Cb_Distance_FuncCAP;
typedef utility::pointer::weak_ptr< Cen_Distance_Func const > Cen_Distance_FuncCAP;
typedef utility::pointer::weak_ptr< CaCbCb_Angle_Func const > CaCbCb_Angle_FuncCAP;
typedef utility::pointer::weak_ptr< NCaCaC_Dihedral_Func const > NCaCaC_Dihedral_FuncCAP;
typedef utility::pointer::weak_ptr< CaCbCbCa_Dihedral_Func const > CaCbCbCa_Dihedral_FuncCAP;

} //disulfides
} //scoring
} //core


#endif //INCLUDED_core_scoring_disulfides_CentroidDisulfidePotential_FWD_HH
