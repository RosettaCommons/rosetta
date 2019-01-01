// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSSingleSet.fwd.hh
/// @brief   forward declaration for PCSSingleSet.hh
/// @details last Modified: 06/22/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_FWD_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <platform/types.hh>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

class PCSSingleSet;

typedef utility::pointer::shared_ptr< PCSSingleSet > PCSSingleSetOP;
typedef utility::pointer::shared_ptr< PCSSingleSet const > PCSSingleSetCOP;
typedef utility::pointer::weak_ptr< PCSSingleSet > PCSSingleSetAP;
typedef utility::pointer::weak_ptr< PCSSingleSet const > PCSSingleSetCAP;

// Forward declaration of friend function
void pcs_erf(platform::Real const *par, int m_dat, void const *data, platform::Real *fvec, int */*info*/);

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSSingleSet_FWD_HH
