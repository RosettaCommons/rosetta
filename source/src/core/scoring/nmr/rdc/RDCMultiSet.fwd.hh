// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/rdc/RDCMultiSet.fwd.hh
/// @brief   forward declaration for RDCMultiSet.hh
/// @details last Modified: 07/27/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_FWD_HH
#define INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <platform/types.hh>

namespace core {
namespace scoring {
namespace nmr {
namespace rdc {

class RDCMultiSet;

typedef utility::pointer::shared_ptr< RDCMultiSet > RDCMultiSetOP;
typedef utility::pointer::shared_ptr< RDCMultiSet const > RDCMultiSetCOP;
typedef utility::pointer::weak_ptr< RDCMultiSet > RDCMultiSetAP;
typedef utility::pointer::weak_ptr< RDCMultiSet const > RDCMultiSetCAP;

// Forward declaration of friend function
void rdc_erf(platform::Real const *par, int m_dat, void const *data, platform::Real *fvec, int */*info*/);

} // namespace rdc
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_rdc_RDCMultiSet_FWD_HH
