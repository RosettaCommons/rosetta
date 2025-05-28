// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/DrugPolishMover.fwd.hh
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_DrugPolishMover_fwd_hh
#define INCLUDED_protocols_drug_design_DrugPolishMover_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace drug_design {

class DrugPolishMover;
typedef utility::pointer::shared_ptr< DrugPolishMover >  DrugPolishMoverOP;
typedef utility::pointer::shared_ptr< DrugPolishMover const >  DrugPolishMoverCOP;


} // namespace drug_design
} // namespace protocols

#endif
