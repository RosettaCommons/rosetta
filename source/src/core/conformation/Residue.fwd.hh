// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/Residue.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Residue_fwd_hh
#define INCLUDED_core_conformation_Residue_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace conformation {

class Residue;

typedef  utility::pointer::weak_ptr< Residue >  ResidueAP;
typedef  utility::pointer::weak_ptr< Residue const >  ResidueCAP;
typedef  utility::pointer::shared_ptr< Residue >  ResidueOP;
typedef  utility::pointer::shared_ptr< Residue const >  ResidueCOP;

typedef  utility::vector1< ResidueOP >  ResidueOPs;
typedef  utility::vector1< ResidueCOP >  ResidueCOPs;
typedef  utility::vector1< ResidueCAP >  ResidueCAPs;

} // namespace conformation
} // namespace core

#endif // INCLUDED_core_conformation_Residue_FWD_HH
