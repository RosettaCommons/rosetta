// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rings/RingConformer.fwd.hh
/// @brief  RingConformer forward declarations header
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rings_RingConformer_fwd_hh
#define INCLUDED_core_chemical_rings_RingConformer_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace rings {

struct RingConformer;

using PoseOP  = utility::pointer::shared_ptr< RingConformer > ;
using PoseCOP = utility::pointer::shared_ptr< RingConformer const >;

} // namespace rings
} // namespace chemical
} // namespace core

#endif // INCLUDED_core_chemical_rings_RingConformer_FWD_HH
