// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/eval/VallFragmentEval.fwd.hh
/// @brief  forward declaration of core::fragment::picking_old::vall::eval::VallFragmentEval
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_fwd_hh
#define INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_fwd_hh


// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


/// @brief  forward declaration of core::fragment::picking_old::vall::eval::VallFragmentEval
class VallFragmentEval;


/// @brief VallFragmentEval owning pointer
typedef utility::pointer::shared_ptr< VallFragmentEval > VallFragmentEvalOP;


/// @brief VallFragmentEval const owning pointer
typedef utility::pointer::shared_ptr< VallFragmentEval const > VallFragmentEvalCOP;


/// @brief VallFragmentEval access pointer
typedef utility::pointer::weak_ptr< VallFragmentEval > VallFragmentEvalAP;


/// @brief VallFragmentEval const access pointer
typedef utility::pointer::weak_ptr< VallFragmentEval const > VallFragmentEvalCAP;


} // eval
} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_eval_VallFragmentEval_FWD_HH */
