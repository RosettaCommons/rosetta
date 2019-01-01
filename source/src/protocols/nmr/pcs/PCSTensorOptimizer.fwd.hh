// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSTensorOptimizer.fwd.hh
/// @brief   forward declaration for PCSTensorOptimizer.hh
/// @details last Modified: 07/05/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_FWD_HH
#define INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace nmr {
namespace pcs {

class PCSTensorOptimizer;

typedef utility::pointer::shared_ptr< PCSTensorOptimizer > PCSTensorOptimizerOP;
typedef utility::pointer::shared_ptr< PCSTensorOptimizer const > PCSTensorOptimizerCOP;
typedef utility::pointer::weak_ptr< PCSTensorOptimizer > PCSTensorOptimizerAP;
typedef utility::pointer::weak_ptr< PCSTensorOptimizer const > PCSTensorOptimizerCAP;


} // namespace pcs
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_FWD_HH
