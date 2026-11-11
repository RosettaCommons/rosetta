// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/PerturbBarrel.fwd.hh
/// @brief  Defines owning pointers for PerturbBarrel mover class.
/// @author Andy Watkins

#ifndef INCLUDED_protocols_beta_barrel_PerturbBarrel_fwd_hh
#define INCLUDED_protocols_beta_barrel_PerturbBarrel_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace beta_barrel {

class PerturbBarrel; // fwd declaration
typedef utility::pointer::shared_ptr< PerturbBarrel > PerturbBarrelOP;
typedef utility::pointer::shared_ptr< PerturbBarrel const > PerturbBarrelCOP;
typedef utility::vector1<PerturbBarrelOP> PerturbBarrelOPs;
typedef utility::vector1<PerturbBarrelCOP> PerturbBarrelCOPs;

} // namespace beta_barrel
} // namespace protocols

#endif // INCLUDED_protocols_beta_barrel_PerturbBarrel_fwd_hh
