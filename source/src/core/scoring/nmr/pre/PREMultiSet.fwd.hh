// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pre/PREMultiSet.fwd.hh
/// @brief   forward declaration for PREMultiSet.hh
/// @details last Modified: 10/12/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pre_PREMultiSet_FWD_HH
#define INCLUDED_core_scoring_nmr_pre_PREMultiSet_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <platform/types.hh>

namespace core {
namespace scoring {
namespace nmr {
namespace pre {

class PREMultiSet;

typedef utility::pointer::shared_ptr< PREMultiSet > PREMultiSetOP;
typedef utility::pointer::shared_ptr< PREMultiSet const > PREMultiSetCOP;
typedef utility::pointer::weak_ptr< PREMultiSet > PREMultiSetAP;
typedef utility::pointer::weak_ptr< PREMultiSet const > PREMultiSetCAP;

// Forward declaration of friend functions
void pre_erf_opt_tau(platform::Real const *par, int m_dat, void const *data, platform::Real *fvec, int */*info*/);
void pre_erf_opt_tau_xyz(platform::Real const *par, int m_dat, void const *data, platform::Real *fvec, int */*info*/);

} // namespace pre
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pre_PREMultiSet_FWD_HH
