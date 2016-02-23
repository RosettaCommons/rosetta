// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////
/// @file PseudocontactShiftData.fwd.hh
///
/// @authorv Christophe Schmitz , Kala Bharath Pilla
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcsTs4_PseudocontactShiftData_fwd_hh
#define INCLUDED_protocols_scoring_methods_pcsTs4_PseudocontactShiftData_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcsTs4 {

class PCS_data_per_lanthanides_Ts4;

class PCS_data_Ts4;

typedef utility::pointer::shared_ptr< PCS_data_Ts4 > PCS_data_Ts4OP;

}//namespace pcsTs4
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
