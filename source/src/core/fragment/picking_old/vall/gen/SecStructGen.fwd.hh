// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/SecStructGen.fwd.hh
/// @brief  forward declaration for core::fragment::picking_old::vall::gen::SecStructGen
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_fwd_hh
#define INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief forward declaration for core::fragment::picking_old::vall::gen::SecStructGen
class SecStructGen;


/// @brief SecStructGen owning pointer
typedef utility::pointer::shared_ptr< SecStructGen > SecStructGenOP;


/// @brief SecStructGen const owning pointer
typedef utility::pointer::shared_ptr< SecStructGen const > SecStructGenCOP;


/// @brief SecStructGen access pointer
typedef utility::pointer::shared_ptr< SecStructGen > SecStructGenAP;


/// @brief SecStructGen const access pointer
typedef utility::pointer::shared_ptr< SecStructGen const > SecStructGenCAP;


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

#endif /* INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_FWD_HH */
