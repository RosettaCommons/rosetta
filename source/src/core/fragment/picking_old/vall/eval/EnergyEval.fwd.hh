// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/EnergyEval.fwd.hh
/// @brief  forward declaration of core::fragment::picking_old::vall::eval::EnergyEval
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_fwd_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_fwd_hh


// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief forward declaration of core::fragment::picking_old::vall::eval::EnergyEval
class EnergyEval;


/// @brief EnergyEval owning pointer
typedef utility::pointer::shared_ptr< EnergyEval > EnergyEvalOP;


/// @brief EnergyEval const owning pointer
typedef utility::pointer::shared_ptr< EnergyEval const > EnergyEvalCOP;


/// @brief EnergyEval access pointer
typedef utility::pointer::weak_ptr< EnergyEval > EnergyEvalAP;


/// @brief EnergyEval const access pointer
typedef utility::pointer::weak_ptr< EnergyEval const > EnergyEvalCAP;


} // eval
} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_EnergyEval_FWD_HH */
