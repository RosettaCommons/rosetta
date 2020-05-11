// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/EPRSpinLabel.fwd.hh
/// @brief  This is a container class specific to the use of double electron-electron resonance data
/// @details This container manages the simulated distance distributions for the deer_decay and
///      deer_distance energy method. It also manages individual electron coordinate ensembles for
///      a given protein, although a reference to that pose is not stored in this object.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_EPRSpinLabel_fwd_hh
#define INCLUDED_core_scoring_epr_deer_EPRSpinLabel_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzVector.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace epr_deer {

class EPRSpinLabel;

typedef utility::pointer::shared_ptr< EPRSpinLabel > EPRSpinLabelOP;
typedef utility::pointer::shared_ptr< EPRSpinLabel const > EPRSpinLabelCOP;

} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
