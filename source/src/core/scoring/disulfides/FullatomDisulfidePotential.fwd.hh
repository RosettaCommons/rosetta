// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/FullatomDisulfidePotential.fwd.hh
/// @brief  Fullatom Disulfide Potential class forward declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_FullatomDisulfidePotential_fwd_hh
#define INCLUDED_core_scoring_disulfides_FullatomDisulfidePotential_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace disulfides {

class FullatomDisulfidePotential;

typedef utility::pointer::shared_ptr< FullatomDisulfidePotential > FullatomDisulfidePotentialOP;
typedef utility::pointer::shared_ptr< FullatomDisulfidePotential const > FullatomDisulfidePotentialCOP;

typedef utility::pointer::weak_ptr< FullatomDisulfidePotential > FullatomDisulfidePotentialAP;
typedef utility::pointer::weak_ptr< FullatomDisulfidePotential const > FullatomDisulfidePotentialCAP;

class CBSG_Dihedral_Func;
class SGSG_Dihedral_Func;
class CB_Angle_Func;
class SG_Dist_Func;

typedef utility::pointer::shared_ptr< CBSG_Dihedral_Func > CBSG_Dihedral_FuncOP;
typedef utility::pointer::shared_ptr< SGSG_Dihedral_Func > SGSG_Dihedral_FuncOP;
typedef utility::pointer::shared_ptr< CB_Angle_Func > CB_Angle_FuncOP;
typedef utility::pointer::shared_ptr< SG_Dist_Func > SG_Dist_FuncOP;


} // namespace disulfides
} // namespace scoring
} // namespace core

#endif
