// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A FilterMover that also calls report() on apply()
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Jul. 30, 2014

#ifndef INCLUDED_protocols_peptide_deriver_FilterReporterMover_fwd_hh
#define INCLUDED_protocols_peptide_deriver_FilterReporterMover_fwd_hh

// Project headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace moves {

class FilterReporterMover;

typedef utility::pointer::shared_ptr<FilterReporterMover> FilterReporterMoverOP;
typedef utility::pointer::shared_ptr<FilterReporterMover const> FilterReporterMoverCOP;

} // namespace peptide_deriver
} // namespace protocols

#endif
