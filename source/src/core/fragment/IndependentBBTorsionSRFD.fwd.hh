// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/IndependentBBTorsionSRFD.fwd.hh
/// @brief  forward declaration for IndependentBBTorsionSRFD
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_IndependentBBTorsionSRFD_fwd_hh
#define INCLUDED_core_fragment_IndependentBBTorsionSRFD_fwd_hh


#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace fragment {


/// @brief forward declaration for IndependentBBTorsionSRFD
class IndependentBBTorsionSRFD;


/// @brief access pointer for IndependentBBTorsionSRFD
typedef utility::pointer::weak_ptr< IndependentBBTorsionSRFD > IndependentBBTorsionSRFD_AP;


/// @brief const access pointer for IndependentBBTorsionSRFD
typedef utility::pointer::weak_ptr< IndependentBBTorsionSRFD const > IndependentBBTorsionSRFD_CAP;


/// @brief owning pointer for IndependentBBTorsionSRFD
typedef utility::pointer::shared_ptr< IndependentBBTorsionSRFD > IndependentBBTorsionSRFD_OP;


/// @brief const owning pointer for IndependentBBTorsionSRFD
typedef utility::pointer::shared_ptr< IndependentBBTorsionSRFD const > IndependentBBTorsionSRFD_COP;


} // namespace fragment
} // namespace core


#endif /* INCLUDED_core_fragment_IndependentBBTorsionSRFD_FWD_HH */
