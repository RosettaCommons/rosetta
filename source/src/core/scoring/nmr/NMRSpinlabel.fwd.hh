// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRSpinlabel.fwd.hh
/// @brief   forward declaration for NMRSpinlabel.hh
/// @details last Modified: 08/31/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRSpinlabel_FWD_HH
#define INCLUDED_core_scoring_nmr_NMRSpinlabel_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace nmr {

class NMRSpinlabel;

typedef utility::pointer::shared_ptr< NMRSpinlabel > NMRSpinlabelOP;
typedef utility::pointer::shared_ptr< NMRSpinlabel const > NMRSpinlabelCOP;
typedef utility::pointer::weak_ptr< NMRSpinlabel > NMRSpinlabelAP;
typedef utility::pointer::weak_ptr< NMRSpinlabel const > NMRSpinlabelCAP;

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRSpinlabel_FWD_HH
